#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
    int myrank, nprocs;
    int data, dataReduzida;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(myrank == 0) data = 5;
    if(myrank != 0) data = 12;

    MPI_Reduce(&data, &dataReduzida, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);   

    printf("Hello from processor %d of %d\n", myrank, nprocs);
    if(myrank == 0) printf("Data reduzida: %d\n", dataReduzida);

    MPI_Finalize();
    return 0;
}