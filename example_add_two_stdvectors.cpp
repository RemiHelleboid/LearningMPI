/*
 *  Prerequisties:
 *     This code runs using an MPI library, either OpenMPI or MPICH2.
 *     These libraries can be installed in either a cluster of computers
 *     or a multicore machine.
 *
 *  How to compile:
 *     mpicc -o vec-add VA-MPI-simple.c
 *
 *  How to execute:
 *     mpirun -np 2 ./vec-add
 *
 *     Note that this executes the code on 2 processes, using the -np command line flag.
 *     See ideas for further exploration of MPI using this code at the end of this file.
 */

#include <stdio.h>   //printf
#include <stdlib.h>  //malloc

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "mpi.h"  // must have a system with an MPI library

/*
 * Definitions
 */
#define MASTER 0    // One process will take care of initialization
#define ARRAY_SIZE  // Size of arrays that will be added together.

/*
 *  In MPI programs, the main function for the program is run on every
 *  process that gets initialized when you start up this code using mpirun.
 */
int main(int argc, char *argv[]) {
    std::size_t number_values = pow(2, 27);
    MPI_Status  status;

    int number_processes;
    int process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

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
            list_1[i] = i;
            list_2[i] = i;
        }
    }

    int                 number_element_per_process = number_values / number_processes;
    std::vector<double> shunk_list_1(number_element_per_process);
    std::vector<double> shunk_list_2(number_element_per_process);
    std::vector<double> shunk_list_sum(number_element_per_process);

    MPI_Scatter(&list_1[0],
                number_element_per_process,
                MPI_DOUBLE,
                &shunk_list_1[0],
                number_element_per_process,
                MPI_DOUBLE,
                MASTER,
                MPI_COMM_WORLD);

    MPI_Scatter(&list_2[0],
                number_element_per_process,
                MPI_DOUBLE,
                &shunk_list_2[0],
                number_element_per_process,
                MPI_DOUBLE,
                MASTER,
                MPI_COMM_WORLD);
    
    if (process_rank == MASTER) {
        list_1.clear();
        list_2.clear();
    }

    for (std::size_t idx_value = 0; idx_value < number_element_per_process; ++idx_value) {
        shunk_list_sum[idx_value] = shunk_list_1[idx_value] * shunk_list_2[idx_value];
    }

    MPI_Gather(&shunk_list_sum[0],
               number_element_per_process,
               MPI_DOUBLE,
               &list_sum[0],
               number_element_per_process,
               MPI_DOUBLE,
               MASTER,
               MPI_COMM_WORLD);

    // if (process_rank == MASTER) {
    //     bool error = false;
    //     for (std::size_t i = 0; i < number_values; i++) {
    //         if (list_sum[i] != i + i) {
    //             error = true;
    //             std::cout << "Error: " << i << " * " << i << " = " << list_sum[i] << std::endl;
    //         }
    //     }
    //     if (error) {
    //         std::cout << "Error: Sum of two vectors is not equal to the sum of the two vectors." << std::endl;
    //     } else {
    //         std::cout << "Success: Sum of two vectors is equal to the sum of the two vectors." << std::endl;
    //     }
    // }

    MPI_Finalize();

    return 0;
}