#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mpi.h"

void poissonSolver(  int numProcs, int myRank,
                    int num_grid_x, int num_grid_y, int num_random );
double f(double x, double y);
void grid( int*, int*,int,int);


const double pi = 3.1415926535897932;
int n1, n2, n3, n4;   // boundary position 
int point = 1;        // where the calculated point is, initializing
                      // 1 means calculation starts from dowm side
int initcoordX[] = {1,99,99,1,50,99,50,1};
int initcoordY[] = {1,99,1,99,1,50,99,50};
int initDirec[] = {1,3,2,4,1,2,3,4};
int sentX[] = {1,-1,0,0,1,0,-1,0};
int sentY[] = {0,0,1,-1,0,1,0,-1};
int lim;
int direcs=1;
        
int main(int argc, char* argv[]) {
    int myRank, numProcs ;
    double startTime, endTime;
    unsigned int semilla;	

    MPI_Init(&argc,&argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
    MPI_Comm_size( MPI_COMM_WORLD, &numProcs );

    // Gather the information needed to parameterize the simulation.

    int num_random = 4000;
    int num_grid_x = 100;
    int num_grid_y = 100;

    semilla=(unsigned int)time((time_t *)NULL);
    MPI_Bcast(&semilla, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    srand((unsigned int)(semilla+131*myRank));

    startTime = MPI_Wtime();
    poissonSolver( numProcs, myRank, num_grid_x, num_grid_y, num_random);
    endTime = MPI_Wtime();

    printf( "The running time (seconds) is %13.8f on processor %d\n",endTime - startTime, myRank );

    MPI_Finalize();
    return 0;
}



// solving the Poisson equation with Dirichlet boundary condition
void poissonSolver( int numProcs, int myRank,
                    int num_grid_x, int num_grid_y, int N ) {

    const int nx = num_grid_x;      // the size of grids
    const int ny = num_grid_y; 
    const double delta = 1. / nx;   // step of x and y
    double ran;
    int ix, iy;


    const int nxx = 100;     // statically allocated size of an array  
    const int nyy = 100;
    double u[nxx+1][nyy+1], x[nxx+1], y[nyy+1];
    int onBoundary[nxx+1][nyy+1]; 

    // cerear u, onBoundary
    int tmpX, tmpY;

    for ( ix=0; ix<=nxx; ix++ ) {
        for ( iy=0; iy<=nyy; iy++ )
            u[ix][iy] = 0.0;
    }

    for ( ix=0; ix<=nx; ix++ ) {
        for ( iy=0; iy<=ny; iy++ )
            onBoundary[ix][iy] = 0;
    }

    // set domain x(0,1) y(0,1)
    for ( ix=0; ix<=nx; ix++ ) x[ix] = (1. * ix) / nx;
    for ( iy=0; iy<=ny; iy++ ) y[iy] = (1. * iy) / ny;

    // set outer boundary condition
    for ( ix=0; ix<=nx; ix++ ) {
        u[ix][0] = 0;    // lower
        onBoundary[ix][0] = 1;
        u[ix][ny] = 0;   // upper
        onBoundary[ix][ny] = 1;
    }

    for ( iy=0; iy<=ny; iy++ ) {
        u[0][iy] = 0;    // left 
        onBoundary[0][iy] = 1;
        u[nx][iy] = 0;   // right
        onBoundary[nx][iy] = 1;
    } 

    // set inner boundary condition
    int inx1 = nx / 4;
    int inx2 = 3 * nx / 4;
    int iny1 = ny / 4;
    int iny2 = 3 * ny / 4;

    for ( ix=inx1; ix<=inx2; ix++ ) {
        u[ix][iny1] = 0;    // lower
        onBoundary[ix][iny1] = 1;
        u[ix][iny2] = 0;   // upper
        onBoundary[ix][iny2] = 1;
    }

    for ( iy=iny1; iy<=iny2; iy++ ) {
        u[inx1][iy] = 0;    // left
        onBoundary[inx1][iy] = 1;
        u[inx2][iy] = 0;   // right
        onBoundary[inx2][iy] = 1;
    }



    // initializing
    ix = initcoordX[myRank];
    iy = initcoordY[myRank];
    point = initDirec[myRank];
    lim = (2*(ny-1)+2*(nx-1))/numProcs;
    n1 = ix; n2=iy;
    n3 = ix+lim*sentX[myRank]-1*sentX[myRank]; n4=iy+lim*sentY[myRank]-1*sentY[myRank];

    if(lim>nx) {
        direcs=2;
        if(n1==1) n3=99,n4=99;
        if(n1==99) n3=1,n4=1;
    }
    
    if(numProcs==1) n1=1,n2=1,n3=99,n4=99;
    //if(myRank==3)printf("%d %d %d %d %d %d\n",n1,n2,n3,n4,ix,iy);
    //return;
    int iix, iiy; 
    int flag = 0;   // (= 1) means calculation done
    double myU;
    while( 1 ) {  

       // calculate u[ix][iy] for a point (ix,iy) 
       double uxy = 0, sqU = 0; 

       if (onBoundary[ix][iy] == 1) {   // u[x][y] evaluated before
            myU = -1000.5;
            flag = 1;
            //printf("FFF: %d %d %d %d %d %d %d %d %d\n",myRank,ix,iy,point,n1,n2,n3,n4);
       }
       else {
            int i;
            for ( i=1; i<=N; i++ ) {   // N random walks
                iix = ix;
                iiy = iy;
                double fxy = f( x[iix], y[iiy] );
                // while loop for one random walk
                while ( 1 ) {
                    ran=(rand())/(double) RAND_MAX;
                    if ( ran >= 0.75 ) iix++;
                    else if ( ran >= 0.50) iiy++; 
                    else if ( ran >= 0.25) iix--;   
                    else iiy--;
                    if ( onBoundary[iix][iiy] == 1 ) break;
                    fxy += f( x[iix], y[iiy] );
                 }
                double tempU = u[iix][iiy] + fxy * delta * delta / 4;
                uxy += tempU;
                sqU += tempU * tempU;
            }
            myU = uxy / N;
            sqU /= N;
            u[ix][iy] = myU;
            onBoundary[ix][iy] = 1;   // set this node as boundary

            // output
            //if ( ix == iy )
            //         printf("ix=%d iy=%d u=%f sigma=%f on proc. %d\n",ix , iy, u[ix][iy], sqrt( (sqU - myU * myU)/N ), myRank);
       }  
         

        MPI_Status status;
        MPI_Request reqs[40];
        /////// PARALELIZACION VA AQUI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // sends u[ix][iy] to the other processors
       int i,nrq;
       nrq=0;
       double data[] = {(double) ix, (double) iy, myU};
       double datar[4][3];
       for ( i = 0; i < numProcs; i++) {
          if (i != myRank) {
             MPI_Isend(data, 3, MPI_DOUBLE, i, 99, MPI_COMM_WORLD,&reqs[nrq]);
             nrq++;
             }
       }
       int srecv=nrq; 
       // receives u[iix][iiy] from the other processors
       for ( i = 0; i < numProcs; i++) {
          if (i != myRank) { 
             MPI_Irecv(datar[i], 3, MPI_DOUBLE, i, 99, MPI_COMM_WORLD,&reqs[nrq]);
             nrq++;
          }  
       }
       
       
       // set ix, iy of next node for the processor
       
        //if(myRank==0)
        //printf("%d %d\n",ix,iy);
       tmpX = ix;
       tmpY = iy;
       grid( &tmpX, &tmpY,myRank,numProcs);
       //if(myRank==4) printf("%d %d %d %d\n",ix,iy,tmpX,tmpY);
       ix = tmpX;
       iy = tmpY;
       
       
       
       int f,index;
       for(i=0;i<nrq;i++){
            f = 0;
            do{
                MPI_Testany(nrq,reqs,&index,&f,&status);
            }
            while(!f);
            if(index>=srecv){
                int s = status.MPI_SOURCE;
                if (datar[s][2] <= -1000) flag = 1;
                 else {   
                    //tmpX = ix;
                    //tmpY = iy;  
                    //grid( &tmpX, &tmpY, i - myRank, 0 );
                    iix = (int)datar[s][0];//tmpX;
                 	iiy = (int)datar[s][1];//tmpY;
                      
                    u[iix][iiy] = datar[s][2];
                    onBoundary[iix][iiy] = 1;
                 }
            }
       }
              
       // check if the whole calculation done, flag = 1
       if (flag == 1) break;

       assert( ((ix < nx+1) && (iy < ny+1))  );
    }
    /* a un archivo de texto */
    int i,j;
    FILE *pFile;
    if (myRank==0)
    {
    pFile = fopen ("resultadoA.txt","w");
    fprintf (pFile, "\t");
    for (i=0; i<=ny; i++) fprintf (pFile, "%.13lf\t",y[i]);
    fprintf (pFile, "\n");
    for (i=0; i<=nx; i++)
        {
        fprintf (pFile, "%.4lf\t",x[i]);
        for (j=0; j<=ny; j++)			
	        fprintf (pFile, " %13.8f\t",u[i][j]);
        fprintf (pFile, "\n");
        }
    fclose(pFile);
    };


    // printf("========  FIN ========= %d\n",myRank);
 

}


// evaluate f(x,y)
double f(double x, double y) {
    return ( 2. * pi * pi * sin(pi*x) * sin(pi*y) /*+ (1-exp(2*x))*/ );
} 


//
void grid(int *x, int *y,int myRank,int nproc) {

     switch ( point ) {
           case 1 :
		    (*x) ++;
		    if ( (*x) >= n3) {
                if(direcs==2 || nproc==1){
                    if((*x)==n3) break;
                    else (*x)--;
                    point=2;
                    (*y)++;
                }else if(nproc>4){
                    if(myRank<4) {n1++;n2++;}
                    else {n3--; n4++; n2++;}
                    (*x) = n1;
                    (*y) = n2;
                }else{
                    n1++; n2++; n4++;n3--;
                    (*x) = n1;
                    (*y) = n2;                
                }
                       
		    }
		    break;
           case 2 :
            (*y)++;
            if ( (*y) >= n4 ) {
               if(nproc==1){
                    if((*y)==n4) break;
                     else (*y)--;
                     point=3;
                      (*x)--; 
               }
               else if(direcs==2){
                point=1;
                n1++; n2++; n3--; n4--;
                (*x) = n1;
                (*y) = n2;
               } else if(nproc>4){
                if(myRank<4) {n1--;n2++;}
                else {n3--; n4--;n1--;}
                (*x) = n1;
                (*y) = n2;
                }else{
                    n1--; n2++; n4--;n3--;
                    (*x) = n1;
                    (*y) = n2;                
                }
            }
            break;
           case 3 :
            (*x) --;
            if(nproc==1 && ((*x)<= n1)){
                //printf("%d %d\n",(*x),n3);
                if((*x)==n1) break;
                else (*x)++;
                point=4;
                (*y)--;        
                 n2++;      
            }
            if(nproc==1) break;
            if ( (*x)<= n3 ) {
                 
                if(direcs==2){
                if((*x)==n3) break;
                else (*x)++;
                point=4;
                (*y)--;
               }else if(nproc>4){
                if(myRank<4) {n1--;n2--;}
                else {n3++; n4--;n2--;}
                (*x) = n1;
                (*y) = n2;
                }else{
                    n1--; n2--; n4--;n3++;
                    (*x) = n1;
                    (*y) = n2;                
                }
	        }       
            break;
           case 4 :
            (*y)--;
             if(nproc==1 && ((*y)< n2)){
                
                point=1;
                (*x)++; (*y)++;
                n1++; n3--; n4--;        
                 //printf("%d %d\n",(*x),n2);   
            }
            if(nproc==1) break;
            if ( (*y) <=  n4 ) {
             if(direcs==2){
                point=3;
                n1--; n2--; n3++; n4++;
                (*x) = n1;
                (*y) = n2;
               }else if(nproc>4){
                if(myRank<4) {n1++;n2--;}
                else {n3++; n4++;n1++;}
                (*x) = n1;
                (*y) = n2;
                }else{
                    n1++; n2--; n4++;n3++;
                    (*x) = n1;
                    (*y) = n2;                
                }               
            }
    }
}
