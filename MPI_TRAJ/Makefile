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

$(OBJ_PATH)basis_r.o: $(SRC_PATH)basis_r.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)awp.o: $(SRC_PATH)awp.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)vmblock.o: $(SRC_PATH)vmblock.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)fgauss.o: $(SRC_PATH)fgauss.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)t_dgls.o: $(SRC_PATH)t_dgls.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH) 

$(OBJ_PATH)psp_pes.o: $(SRC_PATH)psp_pes.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)gear.o: $(SRC_PATH)gear.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)leg_arr.o: $(SRC_PATH)leg_arr.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)matrix.o: $(SRC_PATH)matrix.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)new_matrix.o: $(SRC_PATH)new_matrix.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)co2_ar_dipole.o: $(SRC_PATH)co2_ar_dipole.cpp
	$(CCXX)  $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

objects = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o leg_arr.o co2_ar_dipole.o matrix.o 
OBJ = $(patsubst %,$(OBJ_PATH)%,$(objects))

new_objects = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o leg_arr.o co2_ar_dipole.o new_matrix.o 
NEW_OBJ = $(patsubst %,$(OBJ_PATH)%,$(new_objects))

# Build final binary using MPI compiler wrapper
main: $(OBJ) $(SRC_PATH)main.cpp
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $^ -o $@ $(INC_PATH) 

# using dipole in molecular frame for FFT
new_main: $(NEW_OBJ) $(SRC_PATH)new_main.cpp
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $(TARGET_PATH) $^ -o $@ $(INC_PATH)  
