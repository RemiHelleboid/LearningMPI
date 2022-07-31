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
#include <cmath>
#include <fstream>
#include <iostream>
#include <chrono>
#include <vector>

double compute_eigen_values_random_matrix(double mean, int size) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;
    Eigen::MatrixXd                                matrix1 = Eigen::MatrixXd::Random(size, size);
    Eigen::MatrixXd                                matrix2 = Eigen::MatrixXd::Random(size, size);
    Eigen::MatrixXd                                matrix = matrix1 - matrix2;
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
    std::size_t number_experiments = 10;
    double      min_mean_random    = -5.0;
    double      max_mean_random    = 5.0;
    double      delta_mean         = (max_mean_random - min_mean_random) / (number_values - 1);

    std::cout << "Number of values: " << number_values << std::endl;
    std::cout << "Size of matrix: " << size_matrix << std::endl;
    std::cout << "Delta mean: " << delta_mean << std::endl;
    std::cout << "Min mean: " << min_mean_random << std::endl;
    std::cout << "Max mean: " << max_mean_random << std::endl;

    std::vector<double> list_eigen_values(number_values);

    // Measure time
    auto start = std::chrono::high_resolution_clock::now();
    
    for (std::size_t i = 0; i < number_values; ++i) {
        double sum  = 0.0;
        double mean = min_mean_random + i * delta_mean;
        for (std::size_t j = 0; j < number_experiments; ++j) {
            sum += compute_eigen_values_random_matrix(mean, size_matrix);
        }
        list_eigen_values[i] = sum / number_experiments;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;

    // std::ofstream file_output("output.csv");

    // for (std::size_t i = 0; i < number_values; ++i) {
    //     double mean = min_mean_random + i * delta_mean;
    //     std::cout << "mean: " << mean << " eigen_value: " << list_eigen_values[i] << std::endl;
    //     file_output << mean << " " << list_eigen_values[i] << std::endl;
    // }

    return EXIT_SUCCESS;
}