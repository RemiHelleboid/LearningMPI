#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include <iostream>
#include <cstdlib>
#include <vector>

namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
  mpi::environment env(argc, argv);
  mpi::communicator world;

  std::srand(time(0) + world.rank());
  std::vector<int> all;
  int mine = -1;
  if (world.rank() == 0) {
    all.resize(world.size());
    std::generate(all.begin(), all.end(), std::rand);
  }
  mpi::scatter(world, all, mine, 0);
  for (int r = 0; r < world.size(); ++r) {
    world.barrier();
    if (r == world.rank()) {
      std::cout << "Rank " << r << " got " << mine << '\n';
    }
  }
  return 0;
}
