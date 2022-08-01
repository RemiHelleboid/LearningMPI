/**
 * @file example_hamiltonian.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2022-07-31
 *
 * In this example, the program compute the eigenvalues of random matrix of the form R + a*I.
 * The matrix is actually symetrized, so that we can use the SelfAdjointEigenSolver.
 * Compare to the program example_hamiltonian_mpi.cpp, here we can extract not only the first eigenvalue, but as many as we want.
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

std::vector<double> compute_eigen_values_random_matrix(double mean, int size_matrix, int number_returned_values) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    Eigen::MatrixXd                                matrix1 = Eigen::MatrixXd::Random(size_matrix, size_matrix);
    Eigen::MatrixXd                                matrix2 = Eigen::MatrixXd::Random(size_matrix, size_matrix);
    Eigen::MatrixXd                                matrix  = matrix1 - matrix2;
    matrix                                                 = matrix + mean * Eigen::MatrixXd::Identity(size_matrix, size_matrix);
    // Get the matrix self-adjoint.
    matrix = (matrix + matrix.adjoint()) / 2.0;
    solver.compute(matrix, Eigen::DecompositionOptions::EigenvaluesOnly);
    std::vector<double> eigen_values;
    for (int i = 0; i < number_returned_values; i++) {
        eigen_values.push_back(solver.eigenvalues()(i));
    }
    return eigen_values;
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

    std::size_t size_matrix         = std::atoi(argv[1]);
    std::size_t number_values       = std::atoi(argv[2]);
    std::size_t number_eigen_values = 1;
    double      min_mean_random     = -5.0;
    double      max_mean_random     = 5.0;
    double      delta_mean          = (max_mean_random - min_mean_random) / (number_values - 1);

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
    std::vector<double> list_first_eigenvalues;

    if (process_rank == MASTER) {
        std::cout << "Master process: " << process_rank << std::endl;
        std::cout << "Number of processes: " << number_processes << std::endl;
        std::cout << "Number of values: " << number_values << std::endl;

        list_means.resize(number_values);
        list_first_eigenvalues.resize(number_values * number_eigen_values);

        for (std::size_t i = 0; i < number_values; i++) {
            list_means[i] = min_mean_random + i * delta_mean;
        }
    }

    // Define the number of elements each process will handle.
    int count     = (number_values / number_processes);
    int remainder = (number_values % number_processes);

    std::vector<int> counts_element_per_process(number_processes);
    std::vector<int> displacements_element_per_process(number_processes);
    for (int i = 0; i < number_processes - 1; i++) {
        counts_element_per_process[i]        = count;
        displacements_element_per_process[i] = i * count;
    }
    counts_element_per_process.back()        = (count + remainder);
    displacements_element_per_process.back() = ((number_processes - 1) * count);

    std::vector<double> chunk_list_means;
    chunk_list_means.resize(counts_element_per_process[process_rank]);

    std::cout << "Process " << process_rank << " will handle " << counts_element_per_process[process_rank] << " elements" << std::endl;

    MPI_Scatterv(list_means.data(),
                 counts_element_per_process.data(),
                 displacements_element_per_process.data(),
                 MPI_DOUBLE,
                 chunk_list_means.data(),
                 counts_element_per_process[process_rank],
                 MPI_DOUBLE,
                 MASTER,
                 MPI_COMM_WORLD);

    list_means.clear();

    std::vector<double> chunk_list_first_eigenvalue;
    chunk_list_first_eigenvalue.resize(counts_element_per_process[process_rank] * number_eigen_values);

    for (std::size_t idx_value = 0; idx_value < counts_element_per_process[process_rank]; ++idx_value) {
        std::vector<double> first_eigen_values =
            compute_eigen_values_random_matrix(chunk_list_means[idx_value], size_matrix, number_eigen_values);
        for (std::size_t idx_eigen_value = 0; idx_eigen_value < number_eigen_values; ++idx_eigen_value) {
            chunk_list_first_eigenvalue[idx_value * number_eigen_values + idx_eigen_value] = first_eigen_values[idx_eigen_value];
        }
    }

    std::vector<int> gather_counts_element_per_process(number_processes);
    std::vector<int> gather_displacements_element_per_process(number_processes);
    for (int i = 0; i < number_processes; i++) {
        if (i == process_rank) {
            gather_counts_element_per_process[i] = counts_element_per_process[i] * number_eigen_values;
            gather_displacements_element_per_process[i] *= displacements_element_per_process[i] * number_eigen_values;
        }
    }

    double t3 = MPI_Wtime();
    // MPI_Barrier(MPI_COMM_WORLD);

    double t_gath = MPI_Wtime();
    MPI_Gatherv(chunk_list_first_eigenvalue.data(),
                gather_counts_element_per_process[process_rank],
                MPI_DOUBLE,
                list_first_eigenvalues.data(),
                gather_counts_element_per_process.data(),
                gather_displacements_element_per_process.data(),
                MPI_DOUBLE,
                MASTER,
                MPI_COMM_WORLD);
    chunk_list_first_eigenvalue.clear();
    t_gath = MPI_Wtime() - t_gath;
    std::cout << "Process " << process_rank << " gathered in " << t_gath << " seconds" << std::endl;

    // double t4 = MPI_Wtime();

    // std::cout << "----------------------------------------------------------" << std::endl;
    // std::cout << "process " << process_rank <<  " ----------------------------> total time: " << t4 - t1 << std::endl;
    // std::cout << "----------------------------------------------------------" << std::endl;
    // std::cout << "Export results to file" << std::endl;
    // std::ofstream file_output("output_2.csv");
    // for (std::size_t i = 0; i < number_values; i++) {
    //     file_output << list_means[i] << ',';
    //     for (std::size_t j = 0; j < number_eigen_values; j++) {
    //         file_output << list_first_eigenvalues[i * number_eigen_values + j] << ',';
    //     }
    //     file_output << std::endl;
    // }
    // file_output.close();

    MPI_Finalize();

    return EXIT_SUCCESS;
}