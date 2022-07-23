#include <mpi.h>

#include <iostream>


#define NB_ELEMENTS 128

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank;
    int number_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        int array_1[NB_ELEMENTS];
        int array_2[NB_ELEMENTS];
        for (int i = 0; i < NB_ELEMENTS; i++) {
            array_1[i] = i;
            array_2[i] = -i;
        }


        int result = MPI_Send(&array_1, NB_ELEMENTS, MPI_INT, 1, 0, MPI_COMM_WORLD);
        if (result == MPI_SUCCESS) {
            std::cout << "Rank 0 OK!" << std::endl;
        }
    } else if (rank == 1) {
        int values[NB_ELEMENTS];
        // Receive the value "value", and store it in the buffer value, which is a data of size 1 and of type MPI_INT.
        // The source is rank 0.
        // The tag is 0.
        // The communicator is MPI_COMM_WORLD.        
        int result = MPI_Recv(&values, NB_ELEMENTS, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (result == MPI_SUCCESS && values[17] == 17) {
            std::cout << "Rank 1 OK!" << std::endl;
        }
    }
    MPI_Finalize();
    return 0;
}