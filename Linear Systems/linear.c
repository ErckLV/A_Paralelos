#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

int main(int argc, char *argv[])
{
  float A[10][10]={ {1.01, 0, 0, 0, 0, 0, 0, 0, 0, .1}, 
                    {.1, 1.01, 0, 0, 0, 0, 0, 0, 0, .1},
                    {.1, .1, 1.01, 0, 0, 0, 0, 0, 0, .1},
                    {.1, .1, .1, 1.01, 0, 0, 0, 0, 0, .1},
                    {.1, .1, .1, .1, 1.01, 0, 0, 0, 0, .1},
                    {.1, .1, .1, .1, .1, 1.01, 0, 0, 0, .1},
                    {.1, .1, .1, .1, .1, .1, 1.01, 0, 0, .1},
                    {.1, .1, .1, .1, .1, .1, .1, 1.01, 0, .1},
                    {.1, .1, .1, .1, .1, .1, .1, .1, 1.01, .1},
                    {.1, .1, .1, .1, .1, .1, .1, .1, .1, 1.01} };
  float x[10];
  float b[10]={.1,.1,.1,.1,.1,.1,.1,.1, 1.11, 1.11};                    
  float L[10][10];
  float p[10]={0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  float pp[10][10];
  float v[10][10];
  float x_sum, x_est, r, v_sum, sk;
  int  ult, escog, k;
  int c1,c2,c3, m;
  
  srand((unsigned)time( NULL ));
  // A = I - L
  for (c1=0;c1<10;c1++)   
      for (c2=0;c2<10;c2++)  
         if (c1==c2) L[c1][c2]=1-A[c1][c2]; else L[c1][c2]=-A[c1][c2];
    
         
  for (c1=0;c1<10;c1++) 
      for (c2=0;c2<10;c2++)
         if (fabs(L[c1][c2])>0.0001)
         {  
          //se arma respecto a L
          pp[c1][c2]=(1-p[c1])/(c1<9? c1+2 : 10);
          v[c1][c2]=L[c1][c2]/pp[c1][c2];
         }
         else 
         {
          pp[c1][c2]=0;
          v[c1][c2]=0;
         };
         
  for (c1=0;c1<10;c1++)   
      {for (c2=0;c2<10;c2++)  printf("%6.3f ",pp[c1][c2]);         
       printf("\n");};
         
         
  printf("Cuantos random walk en total? ");
  scanf("%d",&m);
  for (c1=0;c1<10;c1++)
      {
       x_sum=0;
       //iteramos para cada random walk
       for (c2=0;c2<m;c2++)
              {
               v_sum=1;
               ult=c1;
               x_est=b[ult]/p[ult];
               r=(float)rand()/(float)RAND_MAX;
               //genera todo el random walk
               while (r<(1-p[ult]))
                     {
                      sk=0;
                      k=0;
                      //sk camina por filas hasta hallar las probabilidad deseada
                      do
                          {
                           sk=sk+pp[ult][k];
                           k=k+1;
                          }
                      while (r>sk);
                      //escogido antes de absorcion de markov
                      escog=k-1;
                      
                      //printf("%f %f \n",v[ult][escog],v_sum);
                      v_sum=v_sum*v[ult][escog];
                      x_est=v_sum*b[escog]/p[escog];
                      //salta a la sgte fila en la prox iteracion
                      ult=escog;
                      r=(float)rand()/(float)RAND_MAX;
                     };
               x_sum=x_sum+x_est;
              };     
       //aplicacion montecarlo, hallar promedio de todos los x estimados hallados en cada random walk  
       x[c1]=x_sum/m;
      };

  for (c1=0;c1<10;c1++)  printf("x[%d]= %f \n",c1,x[c1]);
  return 0;
}
