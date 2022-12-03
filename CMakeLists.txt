cmake_minimum_required(VERSION 3.15...3.20)

project(
  Armin
  VERSION "0.0.0"
  DESCRIPTION "LearnMPI"
  LANGUAGES C CXX)

include(FetchContent)
include(CheckCXXCompilerFlag)
include(CMakePrintHelpers)

option(ENABLE_BOOST_MPI
       "Build examples using boost::mpi in addition to MPI examples." ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 20)

# Handling build type nicely Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "RelWithDebInfo")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(
    STATUS
      "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "${default_build_type}"
      CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)


include(FetchContent)

# Find Eigen library Find Eigen library
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
  GIT_TAG "3.4.0"
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)

find_package(Eigen3 NO_MODULE)

if(NOT Eigen3_FOUND)
  set(EIGEN_BUILD_DOC OFF)
  set(BUILD_TESTING OFF)
  set(EIGEN_BUILD_PKGCONFIG OFF)
  message("Fetching Eigen3 lib ...")
  add_compile_definitions(EIGEN_USE_LAPACK=0)
  FetchContent_MakeAvailable(Eigen)
endif()

# this variable controls the default value of the options which are normally set
# to ON (useful if you want to configure a minimal version of Gmsh: e.g. "cmake
# -DDEFAULT=0 -DENABLE_POST=1 -DENABLE_PARSER=1")
set(DEFAULT
    ON
    CACHE INTERNAL "Default value for enabled-by-default options")

find_package(MPI)
add_definitions(-DOMPI_SKIP_MPICXX)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_subdirectory(examples)


if(ENABLE_BOOST_MPI)
  # Find boost set(Boost_USE_STATIC_LIBS ON)
  find_package(Boost QUIET REQUIRED COMPONENTS program_options timer mpi
                                               serialization)
  if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost Not found")
  endif()
  add_subdirectory(boostMPI)
endif()

