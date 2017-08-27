cmake_minimum_required (VERSION 2.6)

project(main Fortran)

# source file
#if ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
#    file(GLOB SRC_LIST "*.f90" "gnu/*.f90")
#elseif ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
#    file(GLOB SRC_LIST "*.f90")
#endif()
file(GLOB SRC_LIST "*.f90")


# compiler and compile flag
set (CMAKE_Fortran_FLAGS "-O2 ")
#set (CMAKE_Fortran_FLAGS "-O2 -xAVX -ip ")
#set (CMAKE_Fortran_FLAGS "-O2 -xSSE4.2 -static ")
#set (CMAKE_Fortran_FLAGS "-O2 -xSSE4.2 -unroll ")

# debug
#set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -g -C" )

if ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" )
    set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack -ffree-line-length-none" )
elseif ( "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel" )
    set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl" )
endif()
MESSAGE( "Fortran COMPILER = " ${CMAKE_Fortran_COMPILER_ID} )
MESSAGE( "Fortran COMPILE FLAGS = " ${CMAKE_Fortran_FLAGS} )

# compile
add_executable (main ${SRC_LIST})

