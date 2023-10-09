cmake_minimum_required (VERSION 3.5)
set (CMAKE_CXX_STANDARD 14)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

project (lammps_pot_generator)

include_directories()
 
add_executable (lammps_pot_generator
  printPotential.cpp
  main.cpp
)
target_include_directories(${PROJECT_NAME} PRIVATE
  ${PROJECT_SOURCE_DIR}/
  ${PROJECT_SOURCE_DIR}/../../dependencies/eigen/
)
