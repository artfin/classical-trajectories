# MPI compiler wrapper
MPICC = mpic++ 

# usual compiler
CCXX = g++

CXXFLAGS = -std=c++11 -O3 -lstdc++ -lfftw3 -I /usr/local/include/eigen3 -march=native -mtune=native

# Path for .c, .h and .o Files
SRC_PATH = ./src/
OBJ_PATH = ./obj/
INC_PATH = -I ./include

# Executable name
TARGET = main

all: $(TARGET) test_reader calc_spectra

# Files to compile
objects = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o leg_arr.o matrix.o co2_ar_dipole.o

OBJ = $(patsubst %,$(OBJ_PATH)%,$(objects))
TARGET_PATH = $(patsubst %,$(SRC_PATH)%,$(TARGET).cpp)

# Build .o first
$(OBJ_PATH)%.o: $(SRC_PATH)%.cpp
	@echo [CCXX] $<
	@$(CCXX) $(CXXFLAGS) -o $@ -c $< $(INC_PATH)

# Build final binary using MPI compiler wrapper
$(TARGET): $(OBJ)
	@echo [INFO] Creating Binary Executable [$(TARGET)]
	@$(MPICC) $(CXXFLAGS) $(TARGET_PATH) -o $@ $^ $(INC_PATH) 

test_reader: $(SRC_PATH)test_reader.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

calc_spectra: $(SRC_PATH)calc_spectra.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -r $(OBJ)
