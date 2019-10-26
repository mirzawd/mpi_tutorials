Introduction
===================

This repository contains examples on MPI library in C++ programing language. The purpose of this repository is to provide beginners to advance level training on MPI subroutines.


Dependencies
===================

Mandatory Requirements:
* CMake 2.8.0 (Careful, this will rise soon)
* A C++ compiler (g++ and clang has been tested)
* A Fortran compiler
* A MPI implementation, we described below how to install OpenMPI.

Installing MPI
===================

The simplest way of getting a distribution of MPI is to download it from
the repositories, if this is not possible or the version is too old, you can
get the latest version of OpenMPI from:
     https://www.open-mpi.org/software/ompi/v3.0/
and install it following the steps mention in the documentation.
Very Important, before installing OpenMPI, make sure that the gfortran compiler is installed.


Compiling code
===================

1- Go to the root directory of the source code. 
2- Enter the command  :::: cd build
3- Enter the command  ::::  cmake ..
4- Enter the command  ::::  make -DMPI_BASE_DIR=/path/to/mpi/installation/path/
5- Run the executable :::: mpirun -np x ../bin/src

ENJOY and please do provide your feedback if you have any questions, like to have an improvement in the source code or want to report a bug. 

 
