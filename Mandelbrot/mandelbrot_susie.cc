/**
 *  \file mandelbrot_susie.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

#include <iostream>
#include <cstdlib>
#include <mpi.h>
#include "render.hh"

#define WIDTH 1000
#define HEIGHT 1000

using namespace std;

#define tag  1
#define done 2

//#define debug

int
mandelbrot(double x, double y) {
  int maxit = 511;
  double cx = x;
  double cy = y;
  double newx, newy;

  int it = 0;
  for (it = 0; it < maxit && (x*x + y*y) < 4; ++it) {
    newx = x*x - y*y + cx;
    newy = 2*x*y + cy;
    x = newx;
    y = newy;
  }
  return it;
}

int
main (int argc, char* argv[])
{
	double minX = -2.1;
	double maxX = 0.7;
	double minY = -1.25;
	double maxY = 1.25;

	int height, width;
	if (argc == 3) {
		height = atoi (argv[1]);
		width = atoi (argv[2]);
		assert (height > 0 && width > 0);
	} 
	else 
	{
		fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
		fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
		return -1;
	}

	double it = (maxY - minY)/height;
	double jt = (maxX - minX)/width;
	double x, y;

	int p;
	int rank;
	
	// start MPI
	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &p);		// determine number of processor
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);		// determine processor rank
	
	
	double recvBuff[width * height];
	
	// master
	if (rank == 0)		
	{
		// image handlers
		gil::rgb8_image_t img(height, width);
		auto img_view = gil::view(img);
		
		
		// receive results via gather
		//MPI_GATHER(sendbuf, sendcount, MPI_DOUBLE, recvBuff, width * (height/p), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		//write to image
		int row = 0;
		rank = 1;
		int n = 0;
		
		/*for (rank = 1;rank < p; rank++)
		{
			row = rank + n*p;		
			n++;
			for (int i = 0; i < height; ++i) 
			{
			
				if( (row < height) && (i == row) )
				{
					for (int j = 0; j < width; ++j) 
					{
						img_view(j, i) = render(recvBuff[a + offset]);		
					}
				}
			}
		}
		
		*/
		
		gil::png_write_view("mandelbrot.png", const_view(img));
		
		// tasks all done so send done msg
		for (int proc = 1; proc < p; proc++)
		{
			MPI_Send(0, 0, MPI_DOUBLE, proc, done, MPI_COMM_WORLD);
			#ifdef debug1
				printf("sent done tag to slave %d\n", proc);
			#endif
		} 
	}
	
	// slave
	else 
	{
		slave_proc(minY, minX, width, height, it, jt);
	}
	
	// end of MPI
	MPI_Finalize();
	return 0;	
}

void slave_proc(double minY, double minX, int width, int height, double it, double jt){
	int rank;
	int p;	//number of processors
	double A[width * (height/p)];	// 1D array stores height/p rows that are size width each
	int row = 0;
	int n = 0;
	int y;
	int x;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Status status;

	while(true)
	{
		int offset = 0;
		//calculate result
		y = minY;
		for(rank = 0; rank < p; rank++)
		{
			row = rank + n*p;		//determine row to compute
			n++;
			for (int i = 0; i < height; ++i) {
			
				if( (row < height) && (i == row) )
				{
					x = minX;
					for (int j = 0; j < width; ++j) 
					{
						// calculate row
						A[j + offset] = mandelbrot(x, y)/512.0;
				
						x += jt;
					}
					y += it;
					offset += width;
				}
			}
		}
		
		// receive done tag
		MPI_Recv(&A, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);	
	
		// check for done tag
		if (status.MPI_TAG == done)
		{
			#ifdef debug
			printf("received done tag in slave\n");
			#endif
			
			break;
		}
	}

}

/* eof */
