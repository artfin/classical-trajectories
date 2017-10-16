# MPI compiler wrapper
MPICC = mpic++ 

# usual compiler
CCXX = g++

CXXFLAGS = -std=c++11 -O3 -lstdc++ -march=native -mtune=native -lm

INCLUDE_FFTW = -lfftw3
INCLUDE_EIGEN = -I/usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas

# Path for .c, .h and .o Files
SRC_PATH = ./src/
OBJ_PATH = ./obj/
INC_PATH = -I ./include

# Executable name
TARGET = main

# Files to compile
objects = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o leg_arr.o matrix.o co2_ar_dipole.o # hamiltonian.o 

OBJ = $(patsubst %,$(OBJ_PATH)%,$(objects))
TARGET_PATH = $(patsubst %,$(SRC_PATH)%,$(TARGET).cpp)

# Build .o first
$(OBJ_PATH)%.o: $(SRC_PATH)%.cpp
	$(CCXX) $(CXXFLAGS) $(INC_PATH) $(INCLUDE_EIGEN) -o $@ -c $< 

# Build final binary using MPI compiler wrapper
$(TARGET): $(OBJ)
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $(TARGET_PATH) -o $@.x $^ $(INC_PATH) 

clean:
	rm -r $(OBJ)
