/*  \file mandelbrot_joe.cc
 *
 *  \brief Implement your parallel mandelbrot set in this file.
 */

/* Hikaru Kasai and Jay Patel */


#include <iostream>
#include <cstdlib>
#include "render.hh"
#include "mpi.h"

//#define debug

using namespace std;

#define WIDTH 1000
#define HEIGHT 1000

#define tag  1
#define done 2

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
main(int argc, char* argv[]) {

	double minX = -2.1;
	double maxX = 0.7;
	double minY = -1.25;
	double maxY = 1.25;

	int height, width;
	if (argc == 3) {
		height = atoi (argv[1]);
		width = atoi (argv[2]);
		assert (height > 0 && width > 0);
	} else {
		fprintf (stderr, "usage: %s <height> <width>\n", argv[0]);
		fprintf (stderr, "where <height> and <width> are the dimensions of the image.\n");
		return -1;
  	}
	
	double it = (maxY - minY)/height;
	double jt = (maxX - minX)/width;
	double x, y;
	
	int n;
	int rank;
	double globalArray[height][width]; 	// stores all the data after gathering
	
	
	// start MPI
	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &n);		// determine number of processor
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);		// determine processor rank
	int N = floor(height/n);
	double localArray[N][width]; //local array for storing computed data in each processor
	
	#ifdef debug1
		printf("In rank %d\n", rank);
	#endif
	
	if((rank >= 1) && (rank < n))
	{
		#ifdef debug1
		printf("In rank %d\n", rank);
		printf("number of processors %d\n", rank);
		#endif
		
	 	y = minY + ( (rank - 1) * N * it);		// account for offset
 		for (int row = (rank - 1) * N; row < (rank * N);row++) 
 		{
	    		x = minX + ( (rank - 1) * N * jt);
	    		for (int j = 0; j < width; ++j) 
	    		{
	      			localArray[row][j] = (mandelbrot(x, y)/512.0);		// write pixel to local array
	      			
	      			#ifdef debug1
	      			printf("localarray[%d][%d]: %f\n", row, j, localArray[row][j]);
	      			#endif
	      			
	      			x += jt;
	    		}
	    		y += it;
	  	}
	  	#ifdef debug
	      	printf("finished slave loops========================\n");
	      	#endif
	}
	
	if(rank == 0)
	{
		
		#ifdef debug
			printf("in root==============\n");
		#endif
		
		MPI_Gather(&localArray[0], N*width, MPI_DOUBLE, &globalArray[0], height*width, MPI_DOUBLE,0, MPI_COMM_WORLD);
		#ifdef debug
		printf("finished gather==================\n");
		#endif
		
		// image handlers
		gil::rgb8_image_t img(height, width);
		auto img_view = gil::view(img);
		
		for (int j = 0; j < height; ++j)
		{
			for (int i = 0; i < width; ++i)
			{
				img_view(i,j) = render(globalArray[i][j]);
				#ifdef debug
					printf("Rendered pixel: %f at [%d][%d]\n",globalArray[i][j],i,j);
				#endif	
			}
		}	
		gil::png_write_view("mandelbrot.png", const_view(img));
	}
	MPI_Finalize();
	return 0;
	
} /* End of main() */




