


file(GLOB EXAMPLES *.cpp)
foreach(EXAMPLE ${EXAMPLES})
  get_filename_component(EXAMPLENAME ${EXAMPLE} NAME_WE)
  message("--- File ${EXAMPLENAME}")
  add_executable(${EXAMPLENAME} ${EXAMPLE})
  target_include_directories(${EXAMPLENAME} PUBLIC ${CMAKE_SOURCE_DIR}/src)
  target_link_libraries(${EXAMPLENAME} MPI::MPI_CXX Eigen3::Eigen Boost::serialization)
  target_link_libraries(${EXAMPLENAME} OpenMP::OpenMP_CXX)
endforeach()
