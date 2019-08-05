# 
# Copyright (C) 1996-2016	The SIESTA group
#  This file is distributed under the terms of the
#  GNU General Public License: see COPYING in the top directory
#  or http://www.gnu.org/copyleft/gpl.txt.
# See Docs/Contributors.txt for a list of contributors.
#
#-------------------------------------------------------------------
# arch.make file for gfortran compiler.
# To use this arch.make file you should rename it to
#   arch.make
# or make a sym-link.
# For an explanation of the flags see DOCUMENTED-TEMPLATE.make

.SUFFIXES:
.SUFFIXES: .f .F .o .c .a .f90 .F90

SIESTA_ARCH = unknown

CC = gcc
FPP = $(FC) -E -P -x c
FC = mpif90
#FC_SERIAL = gfortran

FFLAGS = -O2 -fPIC -ftree-vectorize

AR = ar
RANLIB = ranlib

SYS = nag

SP_KIND = 4
DP_KIND = 8
KINDS = $(SP_KIND) $(DP_KIND)



COMP_LIBS = libsiestaLAPACK.a libsiestaBLAS.a libfdict.a
#BLAS_LIBS = -lblas
#LAPACK_LIBS = -llapack
BLACS_LIBS=/usr/lib/libblacsF77init-openmpi.so /usr/lib/libblacsCinit-openmpi.so /usr/lib/libblacs-openmpi.so
SCALAPACK_LIBS=/usr/lib/libscalapack-openmpi.so

#LIB_FLOOK=-L/home/ICN2/aakhtar/Softwares/siesta-4.1-b3/Docs/build/flook/0.7.0/lib 

FPPFLAGS = $(DEFS_PREFIX)-DFC_HAVE_ABORT -DMPI -DFC_HAVE_FLUSH -DFC_HAVE_ABORT -DGRID_DP -DPHI_GRID_SP -DSIESTA__FLOOK

LIBS =  $(NETCDF_LIBS) $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(MPI_LIBS) $(COMP_LIBS) $(BLACS_LIBS) $(LIB_FLOOK) -lflookall -ldl
INCFLAGS += -I/home/ICN2/aakhtar/Softwares/siesta-4.1-b3/Docs/build/flook/0.7.0/include
LDFLAGS +=-L/home/ICN2/aakhtar/Softwares/siesta-4.1-b3/Docs/build/flook/0.7.0/lib -Wl,-rpath=/home/ICN2/aakhtar/Softwares/siesta-4.1-b3/Docs/build/flook/0.7.0/lib


FLOOK_PATH=/path/to/flook/parent

MPI_INTERFACE = libmpi_f90.a
MPI_INCLUDE = .
# Dependency rules ---------

FFLAGS_DEBUG = -g -O1   # your appropriate flags here...

# The atom.f code is very vulnerable. Particularly the Intel compiler
# will make an erroneous compilation of atom.f with high optimization
# levels.
atom.o: atom.F
	$(FC) -c $(FFLAGS_DEBUG) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F) $< 

.c.o:
	$(CC) -c $(CFLAGS) $(INCFLAGS) $(CPPFLAGS) $< 
.F.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_fixed_F)  $< 
.F90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FPPFLAGS) $(FPPFLAGS_free_F90) $< 
.f.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_fixed_f)  $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INCFLAGS) $(FCFLAGS_free_f90)  $<

