 /*  \file mandelbrot_ms.cc
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
	double result[width + 2]; 		// stores calculated pixel 
	double A[width + 2];		// first element is row value and rest are col values 
	
	// start MPI
	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &n);		// determine number of processor
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);		// determine processor rank
	
	
	if (rank == 0)		// master
	{
		// image handlers
		gil::rgb8_image_t img(height, width);
		auto img_view = gil::view(img);
		
		// send n tasks to processor at first
		y = minY;
		int destRank = 1;
		
		for (int row = 0; row < height; ++row) 
		{
			x = minX;
			A[0] = y;	// first element is always row value
			A[1] = row; 	//2nd element is always row number
			for (int col = 0; col < width; ++col)
			{
				A[col + 2] = x;		// +1 offset for first element of row value
				x += jt;
			}
			MPI_Send(A, width + 2, MPI_DOUBLE, destRank, tag, MPI_COMM_WORLD);	// send array contaiing row (y) value and all col(x) values in that row
			destRank++;
			if (destRank >= n)		// only send to n-1 processors at first
			{
				break;
			}
			y += it;
		}
		
		// after initial send, receive and send new task for completed result
		for (int row = n; row < height; ++row) 	
		{
			MPI_Recv(result, width + 2, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
			
			// write each pixel to image
			for (int j = 0; j < width; ++j)
			{
				img_view(j,result[1]) = render(result[j + 2]);	
				#ifdef debug1
				printf("result: %f\n", result[j+2]);
				#endif
			}	
		
			#ifdef debug1
			printf("received in master from slave %d\n", status.MPI_SOURCE);
			#endif	
			
			x = minX;
				
			A[0] = y;	// y starts off from last value in previous send
			A[1] = row;	
			for (int col = 0; col < width; ++col)
			{
				A[col + 2] = x;		// +2 offset for first element of row value and 2nd element of row number
				x += jt;
			}
			MPI_Send(A, width + 2, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD);
			#ifdef debug1
			printf("sent from master to slave %d\n", status.MPI_SOURCE);
			#endif	
		
			y += it;			
		}
		
		// receive remaining tasks from n processors
		for(int proc = 1; proc < n; proc++)
		{
			MPI_Recv(result, width + 2, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
			// write each pixel to image
			for (int j = 0; j < width; ++j)
			{
				img_view(j,result[1]) = render(result[j + 2]);	
				#ifdef debug1
				printf("writing row to image\n");
				#endif
			}	
			
			#ifdef debug1
			printf("received remaining from slave %d\n", status.MPI_SOURCE);
			#endif
				
		}
		gil::png_write_view("mandelbrot.png", const_view(img));
		
		// tasks all done so send done msg
		for (int proc = 1; proc < n; proc++)
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
		slave(width);
	}	
	
	// end of MPI
	MPI_Finalize();
	
	#ifdef debug
	printf("finished finalizing\n");
	#endif
}

void slave(int width)
{
	int rank;
	double result[width + 2]; 		// stores calculated pixel 
	double A[width + 2];	
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	while(true)
	{
		MPI_Recv(A, width + 2, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);	// receive array of size 2
		#ifdef debug1
		printf("received A[0] (row value): %f and A[1] (row #): %f in slave %d\n", A[0], A[1], rank);
		#endif

		// check for done tag
		if (status.MPI_TAG == done)
		{
			#ifdef debug1
			printf("received done tag in slave\n");
			#endif
			
			break;
		}
		else if (status.MPI_TAG == tag)
		{
			double yVal = A[0];		// first element is always y value
			result[0] = A[0];
			result[1] = A[1];
	
			// calculate 1xWidth row
			for (int col = 0; col < width; ++col)
			{
				result[col + 2] = mandelbrot(A[col + 2], yVal)/512.0;	
				#ifdef debug1
				printf("result: %f\n", result[col + 2]);
				#endif	
			}
	
			// send back calculated row of 1xwidth
			MPI_Send(result, width + 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			#ifdef debug1
			printf("sent result[1] %f from slave %d\n", result[1], rank);
			#endif
		}
	}
}
/* eof */
