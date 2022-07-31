// int count = N / size;
// int remainder = N % size;
// int start = rank * count + min(rank, remainder);
// int stop = (rank + 1) * count + min(rank + 1, remainder);

// for (int i = start; i < stop; ++i) { a[i] = DO_SOME_WORK(); }

#include <stdio.h>   //printf
#include <stdlib.h>  //malloc

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "mpi.h"  // must have a system with an MPI library

double complicated_function(double x, double y) {
    double sum = 0.0;
    for (int i = 0; i < 5; ++i) {
        double a     = sqrt(x * x + y * y) * sin(x * y) * cos(x * y) * exp(-x * y);
        double pow_x = pow(x, pow(3.0, 0.75));
        double pow_y = pow(y, pow(3.0, 0.75));
        double b     = pow_x * pow_y;
        sum += a * b;
    }
    return sum;
}

/*
 * Definitions
 */
#define MASTER 0         // One process will take care of initialization
#define ARRAY_SIZE 2048  // Size of arrays that will be added together.

/*
 *  In MPI programs, the main function for the program is run on every
 *  process that gets initialized when you start up this code using mpirun.
 */
int main(int argc, char *argv[]) {
    // std::size_t number_values = pow(2, 4);

    if (argc != 2) {
        std::cout << "Error: missing argument" << std::endl;
        std::cout << "Usage: " << argv[0] << " <number_values>" << std::endl;
        return 1;
    }

    std::size_t number_values = std::atoi(argv[1]);
    MPI_Status  status;

    int number_processes;
    int process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    double t1, t2;
    t1 = MPI_Wtime();

    std::vector<double> list_1;
    std::vector<double> list_2;
    std::vector<double> list_sum;

    if (process_rank == MASTER) {
        std::cout << "Master process: " << process_rank << std::endl;
        std::cout << "Number of processes: " << number_processes << std::endl;
        std::cout << "Number of values: " << number_values << std::endl;

        list_1.resize(number_values);
        list_2.resize(number_values);
        list_sum.resize(number_values);

        for (std::size_t i = 0; i < number_values; i++) {
            list_1[i] = rand() / (double)RAND_MAX * 2 * M_PI;
            list_2[i] = rand() / (double)RAND_MAX * 2 * M_PI;
            ;
        }
    }

    // Define the number of elements each process will handle.
    int count     = number_values / number_processes;
    int remainder = number_values % number_processes;
    int start     = process_rank * count + std::min(process_rank, remainder);
    int stop      = (process_rank + 1) * count + std::min(process_rank + 1, remainder);

    std::vector<int> counts_element_per_process(number_processes);
    std::vector<int> displacements_element_per_process(number_processes);
    for (int i = 0; i < number_processes - 1; i++) {
        counts_element_per_process[i]        = count;
        displacements_element_per_process[i] = i * count;
    }
    counts_element_per_process.back()        = count + remainder;
    displacements_element_per_process.back() = (number_processes - 1) * count;

    std::vector<double> chunk_list_1;
    std::vector<double> chunk_list_2;
    std::vector<double> chunk_list_sum;

    for (int i = 0; i < number_processes; i++) {
        chunk_list_1.resize(counts_element_per_process[i]);
        chunk_list_2.resize(counts_element_per_process[i]);
        chunk_list_sum.resize(counts_element_per_process[i]);
    }

    double time_scatter_start = MPI_Wtime();
    MPI_Scatterv(list_1.data(),
                 counts_element_per_process.data(),
                 displacements_element_per_process.data(),
                 MPI_DOUBLE,
                 chunk_list_1.data(),
                 counts_element_per_process[process_rank],
                 MPI_DOUBLE,
                 MASTER,
                 MPI_COMM_WORLD);

    MPI_Scatterv(list_2.data(),
                 counts_element_per_process.data(),
                 displacements_element_per_process.data(),
                 MPI_DOUBLE,
                 chunk_list_2.data(),
                 counts_element_per_process[process_rank],
                 MPI_DOUBLE,
                 MASTER,
                 MPI_COMM_WORLD);

    double time_scatter_end = MPI_Wtime();
    double time_scatter     = time_scatter_end - time_scatter_start;
    printf("Process %d: scatter time: %f\n", process_rank, time_scatter);

    if (process_rank == MASTER) {
        list_1.clear();
        list_2.clear();
    }

    for (std::size_t idx_value = 0; idx_value < counts_element_per_process[process_rank]; ++idx_value) {
        chunk_list_sum[idx_value] = complicated_function(chunk_list_1[idx_value], chunk_list_2[idx_value]);
    }

    std::cout << "Number of elements handle by process " << process_rank << ": " << counts_element_per_process[process_rank] << std::endl;

    double time_gather_start = MPI_Wtime();
    MPI_Gatherv(chunk_list_sum.data(),
                counts_element_per_process[process_rank],
                MPI_DOUBLE,
                list_sum.data(),
                counts_element_per_process.data(),
                displacements_element_per_process.data(),
                MPI_DOUBLE,
                MASTER,
                MPI_COMM_WORLD);
    double time_gather_end = MPI_Wtime();
    double time_gather     = time_gather_end - time_gather_start;
    printf("Process %d: gather time: %f\n", process_rank, time_gather);
    // if (process_rank == MASTER) {
    //     for (std::size_t i = 0; i < number_values; i++) {
    //         std::cout << list_sum[i] << std::endl;
    //     }
    // }

    t2 = MPI_Wtime();
    printf("Elapsed time is %f\n", t2 - t1);

    MPI_Finalize();

    return 0;
}