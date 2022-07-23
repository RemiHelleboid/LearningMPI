/**
 * @file example_many_matrix_determinants.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief
 * @version 0.1
 * @date 2022-07-23
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

#include "square_matrix.hpp"
#include "utils.hpp"

int main(int argc, char *argv[]) {
    std::cout << "Hello, World!" << std::endl;
    std::cout << "This is a program to compute the determinant of a matrix" << std::endl;
    std::size_t nb_matrices = 20;
    std::vector<double> time_determinant;
    for (std::size_t size_matrix = 8; size_matrix <= 12; ++size_matrix) {
        std::cout << "Size of the matrix: " << size_matrix << std::endl;
        // Start the clock
        auto start = std::chrono::high_resolution_clock::now();

        for (std::size_t i = 0; i < nb_matrices; ++i) {
            std::cout << "\rComputing matrix " << i + 1 << " / " << nb_matrices << std::flush;
            Matrix M(size_matrix);
            M.set_random("normal", 0, 1);
            M.compute_determinant();
        }
        // Stop the clock
        auto                          end     = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "Elapsed time: " << elapsed.count() << " s" << std::endl;
        time_determinant.push_back(elapsed.count());
    }
    std::cout << "Time per sizes: \n" << time_determinant << std::endl; 

    return EXIT_SUCCESS;
}