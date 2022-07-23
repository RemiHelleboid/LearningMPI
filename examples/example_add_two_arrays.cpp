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

#include <iostream>

#include "mpi.h"  // must have a system with an MPI library

/*
 * Definitions
 */
#define MASTER 0         // One process will take care of initialization
#define ARRAY_SIZE 4096  // Size of arrays that will be added together.

/*
 *  In MPI programs, the main function for the program is run on every
 *  process that gets initialized when you start up this code using mpirun.
 */
int main(int argc, char *argv[]) {
    // elements of arrays a and b will be added
    // and placed in array c
    int a[ARRAY_SIZE];
    int b[ARRAY_SIZE];
    int c[ARRAY_SIZE];

    int total_proc;      // total nuber of processes
    int rank;            // rank of each process
    int n_per_proc;      // elements per process
    int n = ARRAY_SIZE;  // number of array elements
    int i;               // loop index

    MPI_Status status;  // not used in this arguably poor example
                        // that is devoid of error checking.

    // 1. Initialization of MPI environment
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &total_proc);
    // 2. Now you know the total number of processes running in parallel
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // 3. Now you know the rank of the current process

    // Smaller arrays that will be held on each separate process

    // 4. We choose process rank 0 to be the root, or master,
    // which will be used to  initialize the full arrays.
    if (rank == MASTER) {
        // initialize arrays a and b with consecutive integer values
        // as a simple example
        for (i = 0; i < n; i++)
            a[i] = i;
        for (i = 0; i < n; i++)
            b[i] = i;
    }

    // All processes take part in the calculations concurrently

    // determine how many elements each process will work on
    n_per_proc = n / total_proc;
    /////// NOTE:
    // In this simple version, the number of processes needs to
    // divide evenly into the number of elements in the array
    ///////////

    // 5. Initialize my smaller subsections of the larger array
    int ap[n_per_proc];
    int bp[n_per_proc];
    int cp[n_per_proc];

    // 6.
    // scattering array a from MASTER node out to the other nodes
    MPI_Scatter(a, n_per_proc, MPI_INT, ap, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD);
    // scattering array b from MASTER node out to the other node
    MPI_Scatter(b, n_per_proc, MPI_INT, bp, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD);

    // 7. Compute the addition of elements in my subsection of the array
    for (i = 0; i < n_per_proc; i++) {
        cp[i] = ap[i] + bp[i];
    }

    // 8. MASTER node gathering array c from the workers
    MPI_Gather(cp, n_per_proc, MPI_INT, c, n_per_proc, MPI_INT, MASTER, MPI_COMM_WORLD);

    // all concurrent processes are finished once they all communicate
    // data back to the master via the gather function.

    // Master process gets to here only when it has been able to gather from all processes
    if (rank == MASTER) {
        // sanity check the result  (a test we would eventually leave out)
        int good = 1;
        for (i = 0; i < n; i++) {
            // printf ("%d ", c[i]);
            if (c[i] != a[i] + b[i]) {
                std::cout << "Wrong addition at index " << i << std::endl;
                good = 0;
                break;
            }
        }
        if (good) {
            std::cout << "All correct!" << std::endl;
        }
    }
    // 9. Terminate MPI Environment and Processes
    MPI_Finalize();

    return 0;
}