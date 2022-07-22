cmake_minimum_required(VERSION 3.15...3.20)

project(
  Armin
  VERSION "0.0.0"
  DESCRIPTION "LearnMPI"
  LANGUAGES C CXX)

include(FetchContent)
include(CheckCXXCompilerFlag)
include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 20)

# this variable controls the default value of the options which are normally set
# to ON (useful if you want to configure a minimal version of Gmsh: e.g. "cmake
# -DDEFAULT=0 -DENABLE_POST=1 -DENABLE_PARSER=1")
set(DEFAULT
    ON
    CACHE INTERNAL "Default value for enabled-by-default options")

find_package(MPI)
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Find boost set(Boost_USE_STATIC_LIBS ON)
find_package(Boost QUIET REQUIRED COMPONENTS program_options timer mpi
                                             serialization)
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost Not found")
endif()


add_executable(example_hello_world example_hello_world.cpp)
target_link_libraries(example_hello_world MPI::MPI_CXX)

add_executable(example_size_rank example_size_rank.cpp)
target_link_libraries(example_size_rank MPI::MPI_CXX Boost::mpi Boost::serialization)


add_executable(example_add_arrays example_send_arrays.cpp)
target_link_libraries(example_add_arrays MPI::MPI_CXX Boost::mpi Boost::serialization)

add_executable(example_broadcast example_brodcast.cpp)
target_link_libraries(example_broadcast MPI::MPI_CXX) 


add_executable(example_add_two_arrays example_add_two_arrays.cpp)
target_link_libraries(example_add_two_arrays MPI::MPI_CXX) 

