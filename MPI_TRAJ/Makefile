# MPI compiler wrapper
MPICC = mpic++ 

# usual compiler
CCXX = g++

CXXFLAGS = -std=c++11 -O3 -lstdc++ -lfftw3 -march=native -mtune=native -lm

INCLUDE_EIGEN = -I/usr/local/include/eigen3

# Path for .c, .h and .o Files
SRC_PATH = ./src/
OBJ_PATH = ./obj/
INC_PATH = -I ./include

# Executable name
TARGET = main

all: $(TARGET) test_reader calc_spectra

# Files to compile
objects = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o leg_arr.o matrix.o co2_ar_dipole.o hamiltonian.o 

OBJ = $(patsubst %,$(OBJ_PATH)%,$(objects))
TARGET_PATH = $(patsubst %,$(SRC_PATH)%,$(TARGET).cpp)

# Build .o first
$(OBJ_PATH)%.o: $(SRC_PATH)%.cpp
	@echo [CCXX] $<
	@$(CCXX) $(CXXFLAGS) $(INC_PATH) $(INCLUDE_EIGEN) -o $@ -c $< 

# Build final binary using MPI compiler wrapper
$(TARGET): $(OBJ)
	@echo [INFO] Creating Binary Executable [$(TARGET)]
	@$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(TARGET_PATH) -o $@ $^ $(INC_PATH) 

test_reader: $(SRC_PATH)test_reader.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

calc_spectra: $(OBJ) $(SRC_PATH)calc_spectra.cpp 
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INC_PATH) $^ -o $@

clean:
	rm -r $(OBJ)
