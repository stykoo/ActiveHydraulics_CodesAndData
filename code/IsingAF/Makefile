# IsingAF, Anti-ferromagnetic Ising model on the triangular lattice.
# Copyright (C) 2023 ENS de Lyon
# Contributor: Alexis Poncet <alexis.poncet@ens-lyon.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

### MKL ###
# Information on how to link MKL can be found on 
# https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html
MKLROOT=/opt/intel/mkl/
MKL_CFLAGS=-DMKL_ILP64 -m64 -I${MKLROOT}/include
MKL_LDFLAGS=${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_sequential.a ${MKLROOT}/lib/libmkl_core.a -lm -ldl
### END MKL ###

CC=g++
CFLAGS=-g -W -Wall -Og -std=c++17 -march=native $(MKL_CFLAGS) 
LDFLAGS=-g -lboost_program_options -pthread -lpthread $(MKL_LDFLAGS)
EXEC=IsingAF
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: main.cpp
	$(CC) -o $@ -c $< $(CFLAGS)

%.o: %.cpp %.h parameters.h
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
		rm -f *.o

mrproper: clean
		rm -rf $(EXEC)
