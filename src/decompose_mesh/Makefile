#=====================================================================
#
#               S p e c f e m 3 D  V e r s i o n  2 . 1
#               ---------------------------------------
#
#          Main authors: Dimitri Komatitsch and Jeroen Tromp
#    Princeton University, USA and University of Pau / CNRS / INRIA
# (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
#                            April 2011
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================
#
# United States Government Sponsorship Acknowledged.
#

# src/decompose_mesh/Makefile.  Generated from Makefile.in by configure.

SHARED = ../shared/

FC = ifort
FCFLAGS = #-g
MPIFC = mpiifort
MPI_INCLUDES = 
MPILIBS = 
FLAGS_CHECK = -O3 -DFORCE_VECTORIZATION -check nobounds -xHost -ftz -assume buffered_io -assume byterecl -align sequence -vec-report0 -std03 -diag-disable 6477 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage
FCFLAGS_f90 = 

FCCOMPILE_CHECK = ${FC} ${FCFLAGS} $(FLAGS_CHECK) -I${SHARED}
MPIFCCOMPILE_CHECK = ${MPIFC} ${FCFLAGS} $(FLAGS_CHECK) -I${SHARED}
FCLINK = $(FCCOMPILE_CHECK)

CC = gcc
CFLAGS = -g -O2 $(CPPFLAGS)
CPPFLAGS =  -I/usr/include/scotch $(COND_MPI_CPPFLAGS)
COND_MPI_CPPFLAGS = -DWITH_MPI
#COND_MPI_CPPFLAGS =

SCOTCH_LIBS = -L"scotch/lib" -lscotch -lscotcherr
### added support for METIS as well, thus uncomment the line below and compile METIS if you want to use it instead of SCOTCH
#SCOTCH_LIBS = -L"scotch/lib" -lscotch -lscotcherr -L./src/decompose_mesh/metis-4.0.3 -L./metis-4.0.3 -lmetis
SCOTCH_INC = -I"scotch/include"

AR = ar
ARFLAGS = cru
RANLIB = ranlib


# E : executables directory
E = ../../bin
# O : objects directory
O = ../../obj/dec

# objects toggled between the parallel and serial version
#COND_MPI_OBJECTS = $O/parallel.o
#COND_MPI_OBJECTS = $O/serial.o

####
#### targets
####

# default targets for the Pyrized version
DEFAULT = xscotch xdecompose_mesh


default: $(DEFAULT)

all: clean default

decompose_mesh: xdecompose_mesh


####
#### rules for executables
####

DEC_OBJECTS = \
	$O/part_decompose_mesh.o \
	$O/decompose_mesh.o \
	$O/fault_scotch.o \
	$O/get_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
        $O/read_value_parameters.shared.o \
        $O/program_decompose_mesh.o

# rules for the pure Fortran version
xdecompose_mesh:: $(DEC_OBJECTS)
	${FCLINK} -o  $E/xdecompose_mesh $(DEC_OBJECTS) $(SCOTCH_INC)  $(SCOTCH_LIBS)



xscotch:
ifeq (1,1)
	echo "Using bundled Scotch"
	$(MAKE) -C scotch/src
else
	echo "Not using bundled Scotch"
endif


clean:
	rm -f $O/*.o *.mod  $E/xdecompose_mesh
####
#### rule to build each .o file below
####

$O/%.shared.o: $(SHARED)%.f90 $(SHARED)constants.h
	${FCCOMPILE_CHECK} -c -o $@ $<

$O/%.cc.o: $(SHARED)%.c ../../config.h
	${CC} -c $(CFLAGS) $(MPI_INCLUDES) -o $@ $(SHARED)$< -I../../

$O/part_decompose_mesh.o: part_decompose_mesh.f90
	${FCCOMPILE_CHECK} -c -o $O/part_decompose_mesh.o part_decompose_mesh.f90 $(SCOTCH_INC)

$O/fault_scotch.o: fault_scotch.f90
	${FCCOMPILE_CHECK} -c -o $O/fault_scotch.o fault_scotch.f90 $(SCOTCH_INC)

$O/decompose_mesh.o: decompose_mesh.F90 part_decompose_mesh.f90 $O/part_decompose_mesh.o $O/fault_scotch.o
	${FCCOMPILE_CHECK} -c -o $O/decompose_mesh.o decompose_mesh.F90 $(SCOTCH_INC)

$O/program_decompose_mesh.o: program_decompose_mesh.f90 $O/part_decompose_mesh.o $O/decompose_mesh.o
	${FCCOMPILE_CHECK} -c -o $O/program_decompose_mesh.o program_decompose_mesh.f90 $(SCOTCH_INC)


