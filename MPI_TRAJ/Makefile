# MPI compiler wrapper
MPICC = mpic++ 

# usual compiler
CCXX = g++

CXXFLAGS = -std=c++11 -O3 -lstdc++ -march=native -mtune=native -lm

INCLUDE_FFTW = -lfftw3
INCLUDE_EIGEN = -I/usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas

######################################################################
INC_PATH_CONSTANTS = -I ./include/CONSTANTS/
######################################################################

######################################################################
# Gear files; ODE solution
#
SRC_PATH_GEAR = ./src/GEAR/
OBJ_PATH_GEAR = ./obj/GEAR/
INC_PATH_GEAR = -I ./include/GEAR/

$(OBJ_PATH_GEAR)basis_r.o: $(SRC_PATH_GEAR)basis_r.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_GEAR)

$(OBJ_PATH_GEAR)awp.o: $(SRC_PATH_GEAR)awp.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_GEAR)

$(OBJ_PATH_GEAR)vmblock.o: $(SRC_PATH_GEAR)vmblock.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_GEAR)

$(OBJ_PATH_GEAR)fgauss.o: $(SRC_PATH_GEAR)fgauss.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_GEAR)

$(OBJ_PATH_GEAR)t_dgls.o: $(SRC_PATH_GEAR)t_dgls.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_GEAR) 

$(OBJ_PATH_GEAR)gear.o: $(SRC_PATH_GEAR)gear.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_GEAR)
######################################################################

######################################################################
# potential energy surfaces
#
SRC_PATH_POT = ./src/POT/
OBJ_PATH_POT = ./obj/POT/
INC_PATH_POT = -I ./include/POT/

$(OBJ_PATH_POT)leg_arr.o: $(SRC_PATH_POT)leg_arr.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_POT)

$(OBJ_PATH_POT)psp_pes.o: $(SRC_PATH_POT)psp_pes.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_POT)
######################################################################

######################################################################
# hamiltonian and it's deirivatives in matrix form
#
SRC_PATH_CO2AR = ./src/CO2AR/
OBJ_PATH_CO2AR = ./obj/CO2AR/
INC_PATH_CO2AR = -I ./include/CO2AR/

$(OBJ_PATH_CO2AR)matrix.o: $(SRC_PATH_CO2AR)matrix.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ $(INC_PATH_CO2AR) $(INC_PATH_DIP) $(INC_PATH_POT)

$(OBJ_PATH_CO2AR)new_matrix.o: $(SRC_PATH_CO2AR)new_matrix.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ $(INC_PATH_CO2AR)
######################################################################

######################################################################
# Dipole moment
#
SRC_PATH_DIP = ./src/DIP/
OBJ_PATH_DIP = ./obj/DIP/
INC_PATH_DIP = -I ./include/DIP/

$(OBJ_PATH_DIP)co2_ar_dipole.o: $(SRC_PATH_DIP)co2_ar_dipole.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_DIP)
######################################################################

######################################################################
# Miscellaneous functions 
#
SRC_PATH_MISC = ./src/MISC/
OBJ_PATH_MISC = ./obj/MISC/
INC_PATH_MISC = -I ./include/MISC/

$(OBJ_PATH_MISC)aux.o: $(SRC_PATH_MISC)aux.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_FFTW) -c $^ -o $@ $(INC_PATH_MISC)

$(OBJ_PATH_MISC)fft.o: $(SRC_PATH_MISC)fft.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_FFTW) -c $^ -o $@ $(INC_PATH_MISC)

$(OBJ_PATH_MISC)gsl_custom.o: $(SRC_PATH_MISC)gsl_custom.cpp 
	$(CCXX) $(CXXFLAGS) $(INCLUDE_GSL) -c $^ -o $@ $(INC_PATH_MISC)
######################################################################

######################################################################
# CO2-Ar trajectory spectrum
#
GEAR_OBJ = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o gear.o 
GEAR_OBJECTS = $(patsubst %,$(OBJ_PATH_GEAR)%,$(GEAR_OBJ))

POT_OBJ = psp_pes.o leg_arr.o 
POT_OBJECTS = $(patsubst %,$(OBJ_PATH_POT)%,$(POT_OBJ))

DIP_OBJ = co2_ar_dipole.o 
DIP_OBJECTS = $(patsubst %,$(OBJ_PATH_DIP)%,$(DIP_OBJ))

CO2AR_OBJ = matrix.o 
CO2AR_OBJECTS = $(patsubst %,$(OBJ_PATH_CO2AR)%,$(CO2AR_OBJ))

MISC_OBJ = aux.o fft.o gsl_custom.o
MISC_OBJECTS = $(patsubst %,$(OBJ_PATH_MISC)%,$(MISC_OBJ))

CO2AR_TRAJ_OBJ = $(GEAR_OBJECTS) $(POT_OBJECTS) $(DIP_OBJECTS) $(CO2AR_OBJECTS) $(MISC_OBJECTS) $(CONSTANTS_OBJECTS)
CO2AR_TRAJ_INC = $(INC_PATH_GEAR) $(INC_PATH_MISC) $(INC_PATH_CO2AR) $(INC_PATH_CONSTANTS)

co2ar_trajs: $(CO2AR_TRAJ_OBJ) $(SRC_PATH_CO2AR)main.cpp
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $^ -o $@ $(CO2AR_TRAJ_INC) 
######################################################################


######################################################################
# TO DO: CO2AR trajectory spectrum with dipole in molecular frame

#new_objects = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o leg_arr.o co2_ar_dipole.o new_matrix.o 
#NEW_OBJ = $(patsubst %,$(OBJ_PATH)%,$(new_objects))

# using dipole in molecular frame for FFT
#new_main: $(NEW_OBJ) $(SRC_PATH)new_main.cpp
	#$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $(TARGET_PATH) $^ -o $@ $(INC_PATH) 
######################################################################



