#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char *argv[])
{
  float A[10][10]={ {1.01, 0, 0, 0, 0, 0, 0, 0, 0, .1}, 
                    {.1, 1.01, 0, 0, 0, 0, 0, 0, 0, .1},
                    {.1, .1, 1.01, 0, 0, 0, 0, 0, 0, .2},
                    {.1, .1, .1, 1.01, 0, 0, 0, 0, 0, .1},
                    {.1, .1, .1, .1, 1.01, 0, 0, 0, 0, .1},
                    {.1, .1, .1, .1, .1, 1.01, 0, 0, 0, .1},
                    {.1, .1, .1, .1, .1, .1, 1.01, 0, 0, .1},
                    {.1, .1, .1, .1, .1, .1, .1, 1.01, 0, .1},
                    {.1, .1, .1, .1, .1, .1, .1, .1, 1.01, .1},
                    {.1, .1, .1, .1, .1, .1, .1, .1, .1, 1.01} };
  float x[10];
  /* la solucion exacta seria {1,-1,1,-1,4,-4,2,-2,1,1}*/
  float b[10]={1.11,-.81,1.21,-.81,4.14,-3.54,2.12,-1.72, 1.11, 1.11};
  float solexacta[10]={1,-1,1,-1,4,-4,2,-2,1,1};
  float L[10][10];
  float p[10]={0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  float pp[10][10];
  float v[10][10];
  double x_sum, x_est, r, v_prod, sk, x_sumtot ,t0, t1;
  long int  ult, escog, k;
  long int c1,c2, m;
  int nprocs, myproc;

  /* MPI initialization */
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
  MPI_Comm_rank( MPI_COMM_WORLD, &myproc );
  
  for (c1=0;c1<10;c1++)   
      for (c2=0;c2<10;c2++)  
         if (c1==c2) L[c1][c2]=1-A[c1][c2]; else L[c1][c2]=-A[c1][c2];
    
         
  for (c1=0;c1<10;c1++) 
      for (c2=0;c2<10;c2++)
         if (fabs(L[c1][c2])>0.00001)
         {  
          pp[c1][c2]=(1-p[c1])/(c1<9? c1+2 : 10);
          v[c1][c2]=L[c1][c2]/pp[c1][c2];
         }
         else 
         {
          pp[c1][c2]=0;
          v[c1][c2]=0;
         };
         
  if( myproc == 0 ) 
    for (c1=0;c1<10;c1++)   
      {for (c2=0;c2<10;c2++)  printf("%6.3f ",v[c1][c2]);         
       printf("\n");};


  srand((unsigned)time( NULL )+991*myproc);         
  if( myproc == 0 )
    { printf("Cuantos random walk en total? ");
      scanf("%ld",&m);
    }

  MPI_Bcast( &m, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Barrier( MPI_COMM_WORLD );
  t0 = MPI_Wtime();

  for (c1=0;c1<10;c1++)
      {
       x_sum=0;
       x_sumtot=0;
       for (c2=0;c2<m;c2+=nprocs)
              {
               v_prod=1;
               ult=c1;
               x_est=b[ult]/p[ult];
               r=(float)rand()/(float)RAND_MAX;
               while (r<(1-p[ult]))
                     {
                      sk=0;
                      k=0;
                      do
                          {
                           sk=sk+pp[ult][k];
                           k=k+1;
                          }
                      while ((r>sk) && (k<11));
                      escog=k-1;
                      
                      v_prod=v_prod*v[ult][escog];
                      x_est=v_prod*b[escog]/p[escog];
                      ult=escog;
                      r=(float)rand()/(float)RAND_MAX;
                     };
               x_sum=x_sum+x_est;
              };
       MPI_Allreduce( &x_sum, &x_sumtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
       x[c1]=x_sumtot/m;
      };

  t1 = MPI_Wtime();
  if( myproc == 0 ) printf("El calculo tomo %f segs en %d nodos\n", t1-t0, nprocs);

  if( myproc == 0 ) for (c1=0;c1<10;c1++)  printf("x[%ld]= %12.8f err=%10.8f\n", c1,x[c1],fabs(x[c1]-solexacta[c1]));

  MPI_Finalize();
  return 0;
}
