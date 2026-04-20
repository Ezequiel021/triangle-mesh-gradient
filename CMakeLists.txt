cmake_minimum_required(VERSION 3.10)
project(gradiente)

find_package(MPI REQUIRED)

add_executable(gradiente paralelo.cpp)
target_link_libraries(gradiente MPI::MPI_CXX)