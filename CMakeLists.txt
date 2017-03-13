cmake_minimum_required (VERSION 2.6)

project(main Fortran)

#set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
#include(FindMKL)
#find_package( MKL REQUIRED )
file(GLOB SOURCE "./*.f90")
set (CMAKE_Fortran_FLAGS "-g -O3 -C -fpp -mkl")

add_executable (main ${SOURCE})

