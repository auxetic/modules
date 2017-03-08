cmake_minimum_required (VERSION 2.6)

project(main Fortran)

file(GLOB SOURCE "./*.f90")
set (CMAKE_Fortran_FLAGS "-g -O3 -C -fpp")

add_executable (main ${SOURCE})

