#include <stdlib.h>
#include <stdio.h>

#include "timer.h"
#include "cuda_utils.h"

typedef float dtype;

#define N_ (8 * 1024 * 1024)
#define MAX_THREADS 256
#define MAX_BLOCKS 64

#define MIN(x,y) ((x < y) ? x : y)

//#define DEBUG


/* return the next power of 2 number that is larger than x */
unsigned int nextPow2( unsigned int x ) {
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}

/* find out # of threads and # thread blocks for a particular kernel */
void getNumBlocksAndThreads(int whichKernel, int n, int maxBlocks, int maxThreads, int &blocks, int &threads)
{
    if (whichKernel < 3)
    {
				/* 1 thread per element */
        threads = (n < maxThreads) ? nextPow2(n) : maxThreads;
        blocks = (n + threads - 1) / threads;
    }
    else
    {
				/* 1 thread per 2 elements */
        threads = (n < maxThreads*2) ? nextPow2((n + 1)/ 2) : maxThreads;
        blocks = (n + (threads * 2 - 1)) / (threads * 2);
    }
		/* limit the total number of threads */
    if (whichKernel == 5)
        blocks = MIN(maxBlocks, blocks);
}

/* special type of reduction to account for floating point error */
dtype reduce_cpu(dtype *data, int n) {
    dtype sum = data[0];
    dtype c = (dtype)0.0;
    for (int i = 1; i < n; i++)
    {
        dtype y = data[i] - c;
        dtype t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}


__global__ void
kernel4(dtype *g_idata, dtype *g_odata, unsigned int n)
{
	__shared__   volatile dtype scratch[MAX_THREADS];
	unsigned int bid = gridDim.x * blockIdx.y + blockIdx.x;
	unsigned int blockDimNew = blockDim.x * 2; 			// since the new blockDim will be halved for the loop iterations
	unsigned int i = bid * blockDimNew + threadIdx.x;

	// do first iteration of sequential for all threads
	if(i < n) {
		scratch[threadIdx.x] = g_idata[i] + g_idata[i + blockDim.x]; 
	} else {
		scratch[threadIdx.x] = 0;
	}
	__syncthreads ();
	
	int warp_size = 32;

	for(unsigned int stride = (blockDim.x/2); stride > warp_size; stride = (stride/2)) {	
		
		if(threadIdx.x < stride) {				// check index range		
			scratch[threadIdx.x] += scratch[threadIdx.x + stride];
		}
		__syncthreads ();
  	}
  	// manually reduce for last warp (unroll to reduce overhead)
  	if(threadIdx.x < warp_size)
  	{
  		if (n > 64){
      			scratch[threadIdx.x] += scratch[threadIdx.x + 32];
      		}
		if (n > 32){
		scratch[threadIdx.x] += scratch[threadIdx.x + 16];
		}
		scratch[threadIdx.x] += scratch[threadIdx.x + 8];
		scratch[threadIdx.x] += scratch[threadIdx.x + 4];
		scratch[threadIdx.x] += scratch[threadIdx.x + 2];
		scratch[threadIdx.x] += scratch[threadIdx.x + 1]; 
  	}

	if(threadIdx.x == 0) {		// copy back to global array
	   g_odata[bid] = scratch[0];
	}
}




int 
main(int argc, char** argv)
{
	int i;

	/* data structure */
	dtype *h_idata, h_odata, h_cpu;
	dtype *d_idata, *d_odata;	

	/* timer */
	struct stopwatch_t* timer = NULL;
	long double t_kernel_4, t_cpu;

	/* which kernel are we running */
	int whichKernel;

	/* number of threads and thread blocks */
	int threads, blocks;

  int N;
  if(argc > 1) {
    N = atoi (argv[1]);
    printf("N: %d\n", N);
  } else {
    N = N_;
    printf("N: %d\n", N);
  }

	/* naive kernel */
	whichKernel = 4;
	getNumBlocksAndThreads (whichKernel, N, MAX_BLOCKS, MAX_THREADS, 
													blocks, threads);

	/* initialize timer */
	stopwatch_init ();
	timer = stopwatch_create ();

	/* allocate memory */
	h_idata = (dtype*) malloc (N * sizeof (dtype));
	CUDA_CHECK_ERROR (cudaMalloc (&d_idata, N * sizeof (dtype)));
	CUDA_CHECK_ERROR (cudaMalloc (&d_odata, blocks * sizeof (dtype)));

	/* Initialize array */
	srand48(time(NULL));
	for(i = 0; i < N; i++) {
		h_idata[i] = drand48() / 100000;
	}
	CUDA_CHECK_ERROR (cudaMemcpy (d_idata, h_idata, N * sizeof (dtype), cudaMemcpyHostToDevice));	//copy to device
	
	/* ================================================== */
	/* GPU kernel */
	dim3 gb(blocks, 1, 1);
	dim3 tb(threads, 1, 1);

	/* warm up */	
	kernel4 <<<gb, tb>>> (d_idata, d_odata, N);
	cudaThreadSynchronize ();

	stopwatch_start (timer);

	/* execute kernel */
	kernel4 <<<gb, tb>>> (d_idata, d_odata, N);
	int s = blocks;
	while(s > 1) {
		threads = 0;
		blocks = 0;
		getNumBlocksAndThreads (whichKernel, s, MAX_BLOCKS, MAX_THREADS, 
														blocks, threads);

		dim3 gb(blocks, 1, 1);
		dim3 tb(threads, 1, 1);

		kernel4 <<<gb, tb>>> (d_odata, d_odata, s);

		s = (s + threads * 2 - 1) / (threads * 2);
	}
	cudaThreadSynchronize ();

	t_kernel_4 = stopwatch_stop (timer);
	fprintf (stdout, "Time to execute unrolled GPU reduction kernel: %Lg secs\n", t_kernel_4);

  double bw = (N * sizeof(dtype)) / (t_kernel_4 * 1e9);
  fprintf (stdout, "Effective bandwidth: %.2lf GB/s\n", bw);

	/* copy result back from GPU */
	CUDA_CHECK_ERROR (cudaMemcpy (&h_odata, d_odata, sizeof (dtype), cudaMemcpyDeviceToHost));
	/* ================================================== */

	/* ================================================== */
	/* CPU kernel */
	stopwatch_start (timer);
	h_cpu = reduce_cpu (h_idata, N);
	t_cpu = stopwatch_stop (timer);
	fprintf (stdout, "Time to execute naive CPU reduction: %Lg secs\n",
         	 t_cpu);
	/* ================================================== */

	if(abs (h_odata - h_cpu) > 1e-5) {
    fprintf(stderr, "FAILURE: GPU: %f  CPU: %f\n", h_odata, h_cpu);
	} else {
    printf("SUCCESS: GPU: %f  CPU: %f\n", h_odata, h_cpu);
	}

	return 0;
}
