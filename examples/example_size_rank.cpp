#include <mpi.h>
#include <stdio.h>
#include <iostream>

int main(int argc, char** argv) {
    int process_Rank, size_Of_Cluster;

    double pi =3.0;
    std::cout << "pi = " << pi << std::endl;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);

    printf("Hello World from process %d of %d\n", process_Rank, size_Of_Cluster);

    MPI_Finalize();

    return 0;
}