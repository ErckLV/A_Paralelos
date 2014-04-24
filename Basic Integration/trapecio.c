#include <stdio.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
   int nprocs, myproc;
   int i, n;
   double pi_25= 3.141592653589793238462643;
   double pi, sumpi, w, x, y, error;
   double t0, t1;

        /* MPI initialization */
   MPI_Init( &argc, &argv );
   MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
   MPI_Comm_rank( MPI_COMM_WORLD, &myproc );

        /* Determine the number of divisions to use */
   if( argc == 2 ) {      /* Take n from the command line argument */
      sscanf( argv[1], "%d", &n );
   } else {   
        if (myproc == 0)
        {
	    printf("Enter the number of intervals: (0 quits) ");fflush(stdout);
	    scanf("%d",&n);	    
        }
   }
   if( myproc == 0 ) printf("Calculating pi using %d divisions\n", n);
        /* Broadcast the number of divisions to all nodes */

   MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD);  /* ??? Is this needed ??? */
        /* Start the timer after a barrier command */
   MPI_Barrier( MPI_COMM_WORLD );
   t0 = MPI_Wtime();

   pi = 0.0;
   sumpi=0.0;
   w = 1.0 / n;
        /* Each processor starts at a different value and computes its
         * own contributions to pi.     */
   for( i=myproc; i<n; i+=nprocs ) {
      x = (double) i / n;
      y = sqrt( 1.0 - x*x );
      if(i!=0 && i!=n-1) y+=y;
      pi += y * w;
   }

   MPI_Allreduce( &pi, &sumpi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   sumpi = sumpi*2;
   error = fabs( sumpi - pi_25 );

   t1 = MPI_Wtime();
   if( myproc == 0 ) {
      printf("The calculated pi = %f (error = %.14f)\n", sumpi, error);
      printf("The calculation took %f seconds on %d nodes\n", t1-t0, nprocs);
   }
   MPI_Finalize();
   return 0;
}
