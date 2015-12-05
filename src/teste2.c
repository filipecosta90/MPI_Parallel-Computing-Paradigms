#include <stdio.h>
#include "mpi.h"

#define NELTS 50000     /* number of elements in the vector */
#define MASTER 0       /* id of the first process */ 
#define FROM_MASTER 1  /* a message type */ 
#define FROM_WORKER 2  /* a message type */ 

MPI_Status status; 

main(int argc, char **argv) {

  double v[NELTS], /* the vector to be summed */
         sum,    /* the overall vector sum */
         ptSum;    /* a partial sum received from a worker */
  int   numProcs, /* number of processes in virtual machine */
        numWorkers, /* number of worker processes */
        rank,   /* our rank */ 
        source,   /* rank of message source */
        dest,   /* rank of message destination */
        mtype,    /* message type */
        aveSz,    /* average partition size */
        extra,    /* number of left over elements */
        offset,   /* offset of partition in v */
        i,      /* simple counter */
        count;    /* the actual number of elements in a partition */


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  numWorkers = numProcs-1; 

  /******* master process ***********/
  if (rank == MASTER) {
//    Initialize the vector
      for (i=0; i<NELTS; i++)
        v[i]= 1;

    // send vector partitions to the workers
    aveSz = NELTS/numWorkers;
    extra = NELTS%numWorkers;
    offset = 0;
    mtype = FROM_MASTER;
    for (dest=1; dest<=numWorkers; dest++) { // start at 1 not 0 (no master)
      count = (dest <= extra) ? aveSz+1 : aveSz;
      MPI_Send(&count,1,MPI_INT,dest,mtype, MPI_COMM_WORLD);
      MPI_Send(&v[offset],count,MPI_DOUBLE,dest,mtype,MPI_COMM_WORLD);
      /***/ printf("Master sent elements %d to %d to rank %d.\n",
          offset,offset+count,dest);
      offset+=count;
    } 


    // wait for results from workers and sum them
    sum=0;
    mtype = FROM_WORKER;
    for (i=1; i<=numWorkers; i++) { // start at 1 not 0 (no master)
      source = i;
      MPI_Recv(&ptSum,1,MPI_DOUBLE,source,mtype,MPI_COMM_WORLD, &status);
      sum+=ptSum;
    } 
    printf("Vector sum is %lf.\n", sum);
  }  // end master


  /************ worker processes *************/
  if (rank > MASTER) {
    // receive out partition from the master
    mtype = FROM_MASTER;
    source = MASTER;
    MPI_Recv(&count,1,MPI_INT,source,mtype,MPI_COMM_WORLD,&status);
    MPI_Recv(&v,count,MPI_DOUBLE,source,mtype,MPI_COMM_WORLD,&status);

    // compute our partial sum
    sum=0;
    for (i=0;i<count;i++)
      sum+=v[i];

    /***/ printf("Slave %d is sending partial sum %lf to master.\n",rank,sum);
    // send our partial sum to the master
    mtype = FROM_WORKER; 
    MPI_Send(&sum,1,MPI_DOUBLE,MASTER,mtype,MPI_COMM_WORLD);
  } // end worker

  MPI_Finalize();
}  // main
