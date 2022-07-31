#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace mpi = boost::mpi;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Error: missing argument" << std::endl;
        std::cout << "Usage: " << argv[0] << " <number_values>" << std::endl;
        return 1;
    }
    std::size_t number_values = std::atoi(argv[1]);

    mpi::environment  env(argc, argv);
    mpi::communicator world;

    std::srand(time(0) + world.rank());
    std::vector<int> list_all_initial_values;
    std::vector<int> list_all_final_values;
    if (world.rank() == 0) {
        list_all_initial_values.resize(number_values);
        for (std::size_t i = 0; i < number_values; i++) {
            list_all_initial_values[i] = rand() % 100;
        }
        list_all_final_values.resize(number_values);
    }
    int count     = number_values / world.size();
    int remainder = number_values % world.size();

    std::vector<int> counts_element_per_process(world.size());
    std::vector<int> displacements_element_per_process(world.size());
    for (int i = 0; i < world.size() - 1; i++) {
        counts_element_per_process[i]        = count;
        displacements_element_per_process[i] = i * count;
    }
    counts_element_per_process.back()        = count + remainder;
    displacements_element_per_process.back() = (world.size() - 1) * count;

    std::cout << "Process " << world.rank() << ": "
              << "nb elements:: " << counts_element_per_process[world.rank()] <<std::endl;

    std::vector<int> chunk_list_1;
    for (int i = 0; i < world.size(); i++) {
        chunk_list_1.resize(counts_element_per_process[i]);
        for (int j = 0; j < counts_element_per_process[i]; j++) {
            chunk_list_1[j] = 2;
        }
    }
    int out_size;
    mpi::scatterv(world,
                  list_all_initial_values,
                  counts_element_per_process,
                  displacements_element_per_process,
                  chunk_list_1.data(),
                  out_size,
                  0);

    for (int i = 0; i < world.size(); i++) {
        world.barrier();
        std::cout << "Process " << world.rank() << ": " << std::endl;
        for (int j = 0; j < chunk_list_1.size(); j++) {
            std::cout << "data " << j << ": " << chunk_list_1[j] << std::endl;
        }
    }

    std::cout << "Process " << world.rank() << " handles " << out_size << " elements" << std::endl;
    std::cout << "Process " << world.rank() << " should handles " << counts_element_per_process[world.rank()] << " elements" << std::endl;
    std::cout << "Process " << world.rank() << " has size " << chunk_list_1.size() << " " << std::endl;

    return 0;
}
