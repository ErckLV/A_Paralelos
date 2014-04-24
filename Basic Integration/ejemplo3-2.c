#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
   int nprocs, myproc;
   int i, n;
   double pi_25= 3.141592653589793238462643;
   double pi, w, x, y, error, acum, sum_acum;
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
	    printf("Enter the number of random points: (0 quits) ");fflush(stdout);
	    scanf("%d",&n);	    
        }

   }
   if( myproc == 0 ) printf("Calculating pi using %d points\n", n);
        /* Broadcast the number of divisions to all nodes */

   MPI_Bcast( &n, 1, MPI_INT, 0, MPI_COMM_WORLD);  /* ??? Is this needed ??? */
        /* Start the timer after a barrier command */
   MPI_Barrier( MPI_COMM_WORLD );
   t0 = MPI_Wtime();

   pi = 0.0;
   sum_acum=0.0;
   acum=0.0;
   srand( (unsigned)time( NULL )+myproc );
   w = 1.0 / n;
        /* Each processor starts at a different value and computes its
         * own contributions to pi.     */
   int n1 = n* (0.75/0.875);
   int n2 = n* (0.125/0.875);
   /**
    Primer cuadrado x = [0, 0,75] y = [0,1]
   */
   for( i=myproc; i<n1; i+=nprocs ) {
      x = (double) (0.75)*rand()/RAND_MAX;
	  y = (double) rand()/RAND_MAX;
	  
        	  
	  if ((x>0.5) && (y<=1.5-x) && (((x*x)+(y*y)) <= 1.0)) acum=acum+1;
	  if ((x>0.5) && (y>1.5-x)) {
	    //Condicion de proyectar el area en el otro pedazo.
	    double nx = 1.5-x, ny = 1.5-y;
	    if( ((nx*nx) + (ny*ny) )<=1.0 ) acum= acum+1;
	  }
      if ( (x<=0.5) && (((x*x)+(y*y)) <= 1.0) ) acum=acum+1;
   }
          
   /**
    Segundo cuadrado x = [0,75, 1] y = [0, 0,5]
    */
    
   for(i=myproc;i<n2;i+=nprocs){
        x = (double) ((0.25)*rand())/RAND_MAX + 0.75;
        y = (double) (0.5)*rand()/RAND_MAX;
        if( ((x*x) + (y*y)) <= 1.0) acum=acum+1;     
   }
   MPI_Allreduce( &acum, &sum_acum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   pi= (((double) sum_acum)*(0.875))/(double) n;
   pi = pi*4;
   error = fabs( pi - pi_25 );

   t1 = MPI_Wtime();
   if( myproc == 0 ) {
      printf("The calculated pi = %f (error = %f)\n", pi, error);
      printf("The calculation took %f seconds on %d nodes\n", t1-t0, nprocs);
   }
   MPI_Finalize();
   return 0;
}
