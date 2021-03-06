#!/bin/bash

# System-specific settings for winisk.math.uwaterloo.ca
# This also doubles for kazan.math.uwaterloo.ca

CC=icc
CXX=icpc
LD=icpc

# System-specific compiler flags
SYSTEM_CFLAGS="-Wall -wd279 -wd383 -wd444 -wd869 -wd981 -wd1418 -wd1419 -wd1572"
SYSTEM_LDFLAGS=

# Compiler flags for debugging
DEBUG_CFLAGS="-g -DBZ_DEBUG"
DEBUG_LDFLAGS=

# Compiler flags for optimization
OPTIM_CFLAGS="-O3 -fp-model fast=2 -ftz"
OPTIM_LDFLAGS=$OPTIM_CFLAGS

# Compiler flags for extra optimization, such as -ip -ipo on icc
EXTRA_OPTIM_CFLAGS="-ip -ipo"
EXTRA_OPTIM_LDFLAGS=$EXTRA_OPTIM_CFLAGS

# Library names/locations/flags for MPI-compilation.  This will
# probably not be necessary on systems with a working mpicc
# alias
MPI_CFLAGS=
MPI_LIB="-lmpi -lmpi++"
MPI_LIBDIR=-L/opt/sgi/mpt/mpt-2.01/lib
MPI_INCDIR=-I/opt/sgi/mpt/mpt-2.01/include

# Library names/locations for LAPACK
LAPACK_LIB="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core"
LAPACK_LIBDIR=
LAPACK_INCDIR=

# Library locations for blitz; leave blank to use system-installed
# or compiled-with-this-package version
BLITZ_LIBDIR=
BLITZ_INCDIR=

# Library locations for fftw
FFTW_LIBDIR=
FFTW_INCDIR=
# Have to disable SSE2 on itanium machines like winisk/kazan
FFTW_OPTIONS=(--disable-sse2 CFLAGS=${OPTIM_CFLAGS[*]})

# Library locations for UMFPACK
UMF_INCDIR=
UMF_LIBDIR=

# Location/library for BLAS
BLAS_LIB=$LAPACK_LIB
BLAS_LIBDIR=
BLAS_INCDIR=

# Boost toolset override
# Works around a boost bug on itanium machines/icc
# see boost ticket #5001
BOOST_TOOLSET=gcc
