/**
 * @file example_serialization_vector_struct.cpp
 * @author remzerrr (remi.helleboid@gmail.com)
 * @brief Scatter a vector of structs.
 * A list of vectors is created by the master process. This list is staggered to the other processes.
 * Each process receives a vector of structs. The norm of each vector is computed and sent back to the master process.
 * @version 0.1
 * @date 2022-08-03
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
    if (argc < 1) {
        std::cout << "Error: missing argument" << std::endl;
        std::cout << "Usage: " << argv[0] << " <number_k_vector>" << std::endl;
        return 1;
    }

    std::size_t number_vector_k = std::atoi(argv[1]);

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

    // Create a vector of structs.
    std::vector<vector_k> list_vector_k;
    if (process_rank == MASTER) {
        list_vector_k.resize(number_vector_k);
        for (int i = 0; i < number_vector_k; ++i) {
            list_vector_k[i].m_kx = i;
            list_vector_k[i].m_ky = -i;
            list_vector_k[i].m_kz = i;
        }
    }

    // Define the number of elements each process will handle.
    int count     = (number_vector_k / number_processes);
    int remainder = (number_vector_k % number_processes);

    std::vector<int> counts_element_per_process(number_processes);
    std::vector<int> displacements_element_per_process(number_processes);
    for (int i = 0; i < number_processes - 1; i++) {
        counts_element_per_process[i]        = count;
        displacements_element_per_process[i] = i * count;
    }
    counts_element_per_process.back()        = (count + remainder);
    displacements_element_per_process.back() = ((number_processes - 1) * count);

    std::vector<vector_k> chunk_vector_of_k;
    chunk_vector_of_k.resize(counts_element_per_process[process_rank]);

    std::cout << "Process " << process_rank << " will handle " << counts_element_per_process[process_rank] << " elements" << std::endl;

    // Scatter the vector of structs.
    MPI_Scatterv(list_vector_k.data(),
                 counts_element_per_process.data(),
                 displacements_element_per_process.data(),
                 k_vector_type,
                 chunk_vector_of_k.data(),
                 counts_element_per_process[process_rank],
                 k_vector_type,
                 MASTER,
                 MPI_COMM_WORLD);

    // Compute the norm of each element.
    std::vector<double> chunk_list_norm_k;
    chunk_list_norm_k.resize(counts_element_per_process[process_rank]);
    for (int i = 0; i < counts_element_per_process[process_rank]; ++i) {
        chunk_list_norm_k[i] = chunk_vector_of_k[i].norm();
    }

    // Gather the vector of structs.
    // We gather the same aount of elements as the number of processes, so we can reuse the same counts and displacement vectors.
    std::vector<double> list_norm_k;
    list_norm_k.resize(number_vector_k);
    MPI_Gatherv(chunk_list_norm_k.data(),
                counts_element_per_process[process_rank],
                MPI_DOUBLE,
                list_norm_k.data(),
                counts_element_per_process.data(),
                displacements_element_per_process.data(),
                MPI_DOUBLE,
                MASTER,
                MPI_COMM_WORLD);

    // Print the vector of structs.
    // std::cout << "Number vector k " << number_vector_k << std::endl;
    // if (process_rank == MASTER) {
    //     for (int i = 0; i < number_vector_k; ++i) {
    //         std::cout << "Vector " << i << ": Norm -> " << list_norm_k[i] << std::endl;
    //     }
    // }


    MPI_Finalize();
    return EXIT_SUCCESS;
}