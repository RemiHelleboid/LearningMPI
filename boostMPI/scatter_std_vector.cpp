#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>


namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
//   mpi::environment env(argc, argv);
//   mpi::communicator world;

//   std::srand(time(0) + world.rank());
//   std::vector<int> all;
//   int mine;
//   if (world.rank() == 0) {
//     all.resize(100);
//     std::generate(all.begin(), all.end(), std::rand);
//     }
//   mpi::scatter(world, all, mine, 5, 0);
//   std::cout << "Rank " << world.rank() << " got " << *(&mine + 1 ) << '\n';

  return 0;
}
