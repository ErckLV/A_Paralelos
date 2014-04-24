#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mpi.h"
#define NMAX 10000

const double PI = 3.1415926535897932;

/* funcion a integrar sobre R2 es func*normalizada del weight */
double func( double x, double y )
{
  return PI *  ( x*x + y*y );
}



/* probability profile of the random numbers , la integral sobre R2 es pi*/
double weight( double x, double y )
{
  return exp( - ( x*x + y*y ) ) ;
}


/* Metropolis algorithm to generate Random       */
/* Numbers in the plane with probability profile */
/* weight(x,y)                                   */
void Metropolis ( double x[], double y[], double (*weight)(double,double),  
                  double delta, int N, int *count )
{
  double xt, yt, xx, yy, eta, ratio;
  int local_count, ran_num_count;
 
  local_count = 0;
  ran_num_count = 0;

  while( ran_num_count < N-1 )
    { 
       local_count++;
       xx = (double)rand()/(double)RAND_MAX;           /* new random numbers 0 to 1 */
       yy = (double)rand()/(double)RAND_MAX;
       xt = x[ran_num_count] + delta * ( 2*xx - 1.0 ); /* in range -1 to 1 */
       yt = y[ran_num_count] + delta * ( 2*yy - 1.0 );
                                                       /* ration of probability profiles */
       ratio = weight( xt, yt ) / weight( x[ran_num_count] , y[ran_num_count]  );   
       eta = (double)rand()/(double)RAND_MAX;          /* random # in 0-1  */
       if ( ratio > eta )                              /* Metropolis selection */
         {
            ran_num_count++;
            x[ran_num_count]  = xt;                    /* selected numbers */
            y[ran_num_count]  = yt;
         }
    }
  *count = *count + local_count;
}


/* demo program to illustrate Monte-Carlo integration using */
/* Gaussian random numbers in the plane */
/*PARALELO*/
int main( int argc, char *argv[] )
{
  int    round, ir, ir_min, count, ran_num, N, N_total, N_input;
  double x[NMAX], y[NMAX], x_local[NMAX], y_local[NMAX], delta;
  double maxx, maxy, minx, miny, tempm;
  int    plot_control;
  double f, sumf, sumf2, local_sumf, local_sumf2,
         total_sumf, total_sumf2, integral, sigma;
  unsigned int semilla;
  int    node, myid, numprocs;
  int    number, nmax;
  MPI_Status recv_status;

                                                /* join the MPI virtual machine */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  
  nmax = NMAX - ( NMAX % ( numprocs-1) );     /* adjust NMAX for # nodes */
  delta = 0.2;                                /* Metropolis delta - size of jump */
  x[0] = 0.5;                                 /* initial position in plane */
  y[0] = 0.6;
  	maxx=x[0];
    minx=x[0];
    maxy=y[0];
    miny=y[0];
	local_sumf=0.0;
	total_sumf=0.0;
	local_sumf2=0.0;
	total_sumf2=0.0;

  if ( myid == 0 )
    {
                                                /* plot control */
/*     plot_control = 0;
     if ( argc == 2 )
       if ( strcmp( argv[1],  "-p") == 0  )*/
         plot_control = 0;

    N_total = 0;                                /* total # of random numbers */
    ran_num = 0;                                /* current number of random numbers */
    count = 1;                                  /* counter for selected numbers */


    f = func(x[0],y[0]);                        /* Monte-Carlo sums */
    sumf = f;
    sumf2 = f*f ;

    semilla=(unsigned int)time((time_t *)NULL);
    srand((unsigned int)(semilla)+131*myid);
                                                /* accuracy improvement */
    fprintf( stderr, "\n Enter # of random number to use ( 0 to quit ): " );
    while ( scanf( "%d", &N_input ) != EOF && N_input > 0 )
    {
                                                /* commensurate with proc # ? */
      if ( N_input % (numprocs-1) != 0 )
             N_input = N_input - ( N_input % (numprocs-1) );
      printf( " # of random numbers used = %d \n", N_input );
                                                /* add to total */
      N_total = N_total + N_input;

      while ( ran_num < N_total )
        {
                                                /* # of random numbers to generate */
          N = N_total - ran_num;
          if ( N  > nmax ) 
             N = nmax;
                                                /* send N to nodes */
          number = N / (numprocs-1) + 1;

                                                /* generate & send random */
                                                /* numbers to each node */
          for ( node=1 ; node<numprocs ; node++ )
            {
              MPI_Ssend(&number, 1, MPI_INT, node, 121, MPI_COMM_WORLD);

              if ( N != 0 )
                {
                                                /* Metropolis Algorithm */
                  Metropolis ( x, y, weight, delta, number , &count );
                                                /* number control */
                  ran_num = ran_num + number - 1;
                                                /* send x & y to nodes */
                                                /* send x & y to nodes */
                  MPI_Ssend( &x[1], number-1, MPI_DOUBLE, node, 123, MPI_COMM_WORLD );
                  MPI_Ssend( &y[1], number-1, MPI_DOUBLE, node, 123, MPI_COMM_WORLD );
                                                /* plot info */
                  if ( plot_control == 1 )
                    {
                       for ( ir=1 ; ir<number ; ir++ )
                         printf( " %f %f\n", x[ir], y[ir] );
                    }
                                                /* ready for next bunch of */
                                                /* random numbers */
                  x[0] = x[number-1];
                  y[0] = y[number-1];
                }    /* if */
            }        /* for */
                                                /* get sums from nodes */
          MPI_Reduce ( &local_sumf, &total_sumf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
          MPI_Reduce ( &local_sumf2, &total_sumf2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
		  tempm=minx;
		  MPI_Allreduce(&tempm,&minx,1,MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		  tempm=miny;
		  MPI_Allreduce(&tempm,&miny,1,MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
		  tempm=maxx;
		  MPI_Allreduce(&tempm,&maxx,1,MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		  tempm=maxy;
		  MPI_Allreduce(&tempm,&maxy,1,MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

                                                /* add to global totals */
                                                /* add to global totals */
          sumf = sumf + total_sumf;
          sumf2 = sumf2 + total_sumf2;

        }   /* while ran_num */

                                                /* finalize Monte-Carlo calc */
        integral = sumf/(double)N_total;
        sigma = sqrt( ( sumf2/(double)N_total - integral*integral ) / (double)N_total );
		fprintf( stderr, " Max x = %f  Min x = %f Max y = %f  Min y = %f\n", maxx,minx,maxy,miny);
        fprintf( stderr, " Integral = %f  +/- sigma = %f \n", integral, sigma );
        fprintf( stderr, "\n Enter # of random number to use ( 0 to quit ): " );

      }  /* while scanf */
                                               /* siganl -> nodes to suicide */
      number = -1;
      for ( node=1 ; node<numprocs ; node++ )
         MPI_Ssend(&number, 1, MPI_INT, node, 121, MPI_COMM_WORLD);

                                                /* Acceptance ratio - should be > ~0.5 */
      fprintf( stderr, "\n Acceptance ratio: %d %d %f\n\n",
                                   N_total, count,  (double)N_total/(double)count );

                                                /* end main code */
      MPI_Finalize();
      exit(1);

    }    /* if my_node */
  else
    {
                                                /* loop over all bunches of */
                                                /* random numbers */
      for ( ; ; )
        {
                                                /* receive signal or # of numbers */
           MPI_Recv( &number, 1, MPI_INT, 0, 121,
                               MPI_COMM_WORLD, &recv_status);
                                                /* end code - suicide */
           if ( number <0 )
             {
                MPI_Finalize();
                exit(0);
             }
           else if ( number > 0 )              /* perform MC integral */
             {
               MPI_Recv( &x_local, number-1, MPI_DOUBLE, 0, 123,
                               MPI_COMM_WORLD, &recv_status);
               MPI_Recv( &y_local, number-1, MPI_DOUBLE, 0, 123,
                               MPI_COMM_WORLD, &recv_status);
                                                /* Monte-Carlo integral */
               local_sumf = 0.0;
               local_sumf2 = 0.0;
               for ( ir=0 ; ir<number-1 ; ir++ )
                 {
                   f =  func( x_local[ir], y_local[ir] );
                   local_sumf = local_sumf + f;
                   local_sumf2 = local_sumf2 + f*f;
				   if (x_local[ir]<minx) minx=x_local[ir];
				   if (x_local[ir]>maxx) maxx=x_local[ir];
			       if (y_local[ir]<miny) miny=y_local[ir];
			       if (y_local[ir]>maxy) maxy=y_local[ir];
                 }
                                                /* send partial sums to node 0 */
               MPI_Reduce ( &local_sumf, &total_sumf, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
               MPI_Reduce ( &local_sumf2, &total_sumf2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
			   tempm=minx;
			   MPI_Allreduce(&tempm,&minx,1,MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
			   tempm=miny;
			   MPI_Allreduce(&tempm,&miny,1,MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
			   tempm=maxx;
			   MPI_Allreduce(&tempm,&maxx,1,MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			   tempm=maxy;
			   MPI_Allreduce(&tempm,&maxy,1,MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
             }
        }

    }
                                                /* Acceptance ratio - should be > ~0.5 */
  if ( myid == 0 )
     fprintf( stderr, "\n Acceptance ratio: %f\n\n", 
           (double)N_total/(double)count );

                                                /* end of code */
  MPI_Finalize();
  return 0;
}
