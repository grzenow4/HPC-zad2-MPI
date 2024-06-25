cmake_minimum_required(VERSION 3.10)
SET(CMAKE_CXX_COMPILER mpicxx)

project(HPC-zad2-MPI CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-g -O3 -Wall -Wextra -DOMPI_SKIP_MPICXX")

add_executable(matmul matmul.cpp grid.cpp matrix.cpp utils.cpp)
