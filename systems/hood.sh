#!/bin/bash

# System-specific settings for hood.math

CC=icc
CXX=icpc
LD=icpc

# System-specific compiler flags
SYSTEM_CFLAGS=

SYSTEM_LDFLAGS=

# Compiler flags for debugging
DEBUG_CFLAGS="-g -DBZ_DEBUG"
DEBUG_LDFLAGS=

# Compiler flags for optimization
OPTIM_CFLAGS="-O3 -fp-model fast=2"
OPTIM_LDFLAGS=$OPTIM_CFLAGS

# Compiler flags for extra optimization, such as -ip -ipo on icc
EXTRA_OPTIM_CFLAGS="-ip -ipo"
EXTRA_OPTIM_LDFLAGS=$EXTRA_OPTIM_CFLAGS

# Library names/locations/flags for MPI-compilation.  This will
# probably not be necessary on systems with a working mpicc
# alias
MPICXX=icpc
MPI_CFLAGS=
MPI_LIBDIR="-L${MPI_LIB}"
MPI_LIB="-lmpi"
MPI_INCDIR="-I${MPI_INCLUDE}"

# Library names/locations for LAPACK
LAPACK_LIB="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread"
LAPACK_LIBDIR=
LAPACK_INCDIR=

# Library locations for blitz; leave blank to use system-installed
# or compiled-with-this-package version
BLITZ_LIBDIR=
BLITZ_INCDIR=

# Library locations for fftw
FFTW_LIBDIR=
FFTW_INCDIR=

# Library locations for UMFPACK
UMF_INCDIR=
UMF_LIBDIR=

# Location/library for BLAS
# The sharcnet clusters are strange, and their compiler
# script includes the blas libraries with -llapack
BLAS_LIB=$LAPACK_LIB
BLAS_LIBDIR=$LAPACK_LIBDIR
BLAS_INCDIR=

