#include<stdio.h>
#include<string.h>
#include "mpi.h"
#define BUFSIZE 128
#define TAG 0

int main(int argc, char * argv[]){
    printf("%d\n",argc);
    printf("%s\n",argv[0]);
    char idstr[32],buff[BUFSIZE];
    int numprocs,myid,i;
    MPI_Status stat;
    MPI_Init(&argc, &argv); //empieza la paralelizacion
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); //mpi_comm_worl -> devuelve cantidad de procesadores actuales
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); //devuelve el rank del procesador (un id)
    
    //printf("Hola soy el core nro %d de %d\n",myid,numprocs);
    if(myid==0){
        printf("%d hay %d procesadores\n",myid,numprocs);
        for(i=1;i<numprocs;i++){
            sprintf(buff,"Hello %d!",i);
            MPI_Send(buff,BUFSIZE,MPI_CHAR,i,TAG,MPI_COMM_WORLD);
        }
        for(i=1;i<numprocs;i++){
            MPI_Recv(buff,BUFSIZE,MPI_CHAR,i,TAG,MPI_COMM_WORLD,&stat);
            printf("%d :  %s!\n",myid,buff);
        }
    }else{
        MPI_Recv(buff,BUFSIZE,MPI_CHAR,0,TAG,MPI_COMM_WORLD,&stat);
        sprintf(idstr, "Processor %d ", myid);
        strncat(buff, idstr, BUFSIZE-1);
        strncat(buff, "reporting for duty\n", BUFSIZE-1);
        MPI_Send(buff, BUFSIZE, MPI_CHAR, 0, TAG, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}
