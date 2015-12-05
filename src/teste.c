#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
  int rank, numprocs;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  printf("Hello World from proces %d of %d.\n",rank,numprocs);
  MPI_Finalize();
}
