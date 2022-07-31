/**
 * @file example_hamiltonian.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-31
 *
 * In this example, the program is run on every process that gets initialized when you start up this code using mpirun.
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <mpi.h>     // must have a system with an MPI library
#include <stdio.h>   //printf
#include <stdlib.h>  //malloc

#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#define MASTER 0  // One process will take care of initialization

double compute_eigen_values_random_matrix(double mean, int size) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    Eigen::MatrixXd                                matrix1 = Eigen::MatrixXd::Random(size, size);
    Eigen::MatrixXd                                matrix2 = Eigen::MatrixXd::Random(size, size);
    Eigen::MatrixXd                                matrix  = matrix1 - matrix2;
    matrix                                                 = matrix + mean * Eigen::MatrixXd::Identity(size, size);
    // Get the matrix self-adjoint.
    matrix = (matrix + matrix.adjoint()) / 2.0;
    // std::cout << "Matrix: " << std::endl;
    // std::cout << matrix << std::endl;
    solver.compute(matrix, Eigen::DecompositionOptions::EigenvaluesOnly);
    return solver.eigenvalues()(0);
}

/**
 * @brief main function with serial computation.
 *
 * @param argc
 * @param argv
 * @return int
 */
int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cout << "Error: missing argument" << std::endl;
        std::cout << "Usage: " << argv[0] << " <mean> <size>" << std::endl;
        return 1;
    }

    std::size_t size_matrix        = std::atoi(argv[1]);
    std::size_t number_values      = std::atoi(argv[2]);
    double      min_mean_random    = -5.0;
    double      max_mean_random    = 5.0;
    double      delta_mean         = (max_mean_random - min_mean_random) / (number_values - 1);

    // Initialize the MPI environment.
    MPI_Status status;

    int number_processes;
    int process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
    double t1, t2;
    t1 = MPI_Wtime();

    std::vector<double> list_means;
    std::vector<double> list_first_eigenvalue;

    if (process_rank == MASTER) {
        std::cout << "Master process: " << process_rank << std::endl;
        std::cout << "Number of processes: " << number_processes << std::endl;
        std::cout << "Number of values: " << number_values << std::endl;

        list_means.resize(number_values);
        list_first_eigenvalue.resize(number_values);

        for (std::size_t i = 0; i < number_values; i++) {
            list_means[i] = min_mean_random + i * delta_mean;
        }
    }

    // Define the number of elements each process will handle.
    int count     = number_values / number_processes;
    int remainder = number_values % number_processes;

    std::vector<int> counts_element_per_process(number_processes);
    std::vector<int> displacements_element_per_process(number_processes);
    for (int i = 0; i < number_processes - 1; i++) {
        counts_element_per_process[i]        = count;
        displacements_element_per_process[i] = i * count;
    }
    counts_element_per_process.back()        = count + remainder;
    displacements_element_per_process.back() = (number_processes - 1) * count;

    std::vector<double> chunk_list_means;
    for (int i = 0; i < number_processes; i++) {
        chunk_list_means.resize(counts_element_per_process[i]);
    }
    MPI_Scatterv(list_means.data(), counts_element_per_process.data(), displacements_element_per_process.data(),
                 MPI_DOUBLE, chunk_list_means.data(), counts_element_per_process[process_rank], MPI_DOUBLE, MASTER,
                 MPI_COMM_WORLD);

    std::vector<double> chunk_list_first_eigenvalue;
    for (int i = 0; i < number_processes; i++) {
        chunk_list_first_eigenvalue.resize(counts_element_per_process[i]);
    }

    std::cout << "Process " << process_rank << ": " << chunk_list_means.size() << std::endl;

    for (std::size_t idx_value = 0; idx_value < counts_element_per_process[process_rank]; ++idx_value) {
        chunk_list_first_eigenvalue[idx_value] = compute_eigen_values_random_matrix(chunk_list_means[idx_value],
                                                                                    size_matrix);
        // std::cout << "Process " << process_rank << ": " << chunk_list_first_eigenvalue[idx_value] << std::endl;
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    double t3 = MPI_Wtime();
    MPI_Gatherv(chunk_list_first_eigenvalue.data(), counts_element_per_process[process_rank], MPI_DOUBLE,
                list_first_eigenvalue.data(), counts_element_per_process.data(), displacements_element_per_process.data(),
                MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    double t4 = MPI_Wtime();
    // std::cout << "Process gather time " << process_rank << ": " << t4 - t3 << std::endl;
    // if (process_rank == MASTER) {
    //     std::cout << "List of first eigenvalues: " << std::endl;
    //     for (std::size_t i = 0; i < number_values; i++) {
    //         std::cout << "Mean: " << list_means[i] << " First eigenvalue: " << list_first_eigenvalue[i] << std::endl;
    //     }
    // }

    if (process_rank == MASTER) {
        std::cout << "----------------------------------------------------------" << std::endl;
        std::cout << "total time: " << t4 - t1 << std::endl; 
        std::cout << "----------------------------------------------------------" << std::endl;
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}