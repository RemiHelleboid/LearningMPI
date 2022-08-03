/**
 * @file example_serialization.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief Example of serialization of a struct that contains a vector of double.
 * @version 0.1
 * @date 2022-08-02
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <mpi.h>     // must have a system with an MPI library
#include <stdio.h>   //printf
#include <stdlib.h>  //malloc

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#define MASTER 0  // One process will take care of initialization

// Define the structure to be serialized. 
// It is a simple 3D vector of double. 
// We define a unique method, which computes the norm of the vector.
typedef struct vector_k {
    double m_kx;
    double m_ky;
    double m_kz;

    double norm() const { return std::sqrt(m_kx * m_kx + m_ky * m_ky + m_kz * m_kz); }
} vector_k;


int main(int argc, char** argv) {
    // Initialize the MPI environment.
    MPI_Status status;

    int number_processes;
    int process_rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    // Create a new MPI type for the struct k_vector.
    MPI_Datatype k_vector_type;
    const int    number_item_k_vector = 3;
    MPI_Datatype type[3]              = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    int          block_lengths[3]     = {1, 1, 1};

    MPI_Aint offsets[3];
    offsets[0] = offsetof(vector_k, m_kx);
    offsets[1] = offsetof(vector_k, m_ky);
    offsets[2] = offsetof(vector_k, m_kz);

    MPI_Type_create_struct(number_item_k_vector, block_lengths, offsets, type, &k_vector_type);
    MPI_Type_commit(&k_vector_type);
    // END Creating a new MPI type for the struct k_vector.

    // Creating a vector_k instant to do a small test.
    vector_k unitary_vector;
    if (process_rank == MASTER) {
        // Initialize the vector_k only in the master process.
        unitary_vector.m_kx = 1.0 / 2.0;
        unitary_vector.m_ky = 1.0 / 2.0;
        unitary_vector.m_kz = 1.0 / 2.0;
    }

    // Broadcast the vector_k to all processes.
    MPI_Bcast(&unitary_vector, 1, k_vector_type, MASTER, MPI_COMM_WORLD);

    // Modify the vector_k in all the processes.
    unitary_vector.m_kx *= process_rank;
    unitary_vector.m_ky *= process_rank;
    unitary_vector.m_kz *= process_rank;

    // Show the result.
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Process " << process_rank << ": " << unitary_vector.m_kx << " " << unitary_vector.m_ky << " " << unitary_vector.m_kz
              << std::endl;
    
    // Print the norm of the vector_k.
    std::cout << "Process " << process_rank << ": norm of the vector: " << unitary_vector.norm() << std::endl;

    return EXIT_SUCCESS;
}