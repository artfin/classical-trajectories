# MPI compiler wrapper
MPICC = mpic++ 

# usual compiler
CCXX = g++

CXXFLAGS_NOOPT = -std=c++11 
CXXFLAGS = -std=c++11 -O3 -lstdc++ -lm

INCLUDE_FFTW = -lfftw3
INCLUDE_EIGEN = -I/usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas
INCLUDE_GNUPLOT_IO = -I /home/artfin/Desktop/repos/gnuplot-iostream  

OPENGL_LFLAGS = -lglut -lGL -lGLEW -lGLU

LINK_BOOST = -lboost_iostreams -lboost_system -lm

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

# 1 additional parameter
$(OBJ_PATH_GEAR)gear.o: $(SRC_PATH_GEAR)gear.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_GEAR)

# no additional parameters
$(OBJ_PATH_GEAR)gear_diatom.o: $(SRC_PATH_GEAR)gear_diatom.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_GEAR)
######################################################################

######################################################################
# path to gear objects
GEAR_OBJ = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o gear.o 
GEAR_OBJECTS = $(patsubst %,$(OBJ_PATH_GEAR)%,$(GEAR_OBJ))

GEAR_OBJ_DIATOM = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o gear_diatom.o
GEAR_OBJECTS_DIATOM = $(patsubst %,$(OBJ_PATH_GEAR)%,$(GEAR_OBJ_DIATOM))
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

$(OBJ_PATH_POT)ar_he_pes.o: $(SRC_PATH_POT)ar_he_pes.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_POT)

$(OBJ_PATH_POT)ar_he_pes_derivative.o: $(SRC_PATH_POT)ar_he_pes_derivative.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_POT)
######################################################################

######################################################################
# hamiltonian and it's derivatives for co2-ar 
#
SRC_PATH_CO2AR = ./src/CO2AR/
OBJ_PATH_CO2AR = ./obj/CO2AR/
INC_PATH_CO2AR = -I ./include/CO2AR/

$(OBJ_PATH_CO2AR)matrix.o: $(SRC_PATH_CO2AR)matrix.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ $(INC_PATH_CO2AR) $(INC_PATH_DIP) $(INC_PATH_POT)

$(OBJ_PATH_CO2AR)new_matrix.o: $(SRC_PATH_CO2AR)new_matrix.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ $(INC_PATH_CO2AR)

$(OBJ_PATH_CO2AR)matrix_euler.o: $(SRC_PATH_CO2AR)matrix_euler.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ $(INC_PATH_CO2AR) $(INC_PATH_POT) $(INC_PATH_DIP)
######################################################################

######################################################################
# hamiltonian and it's derivatives for diatomics
#
SRC_PATH_DIATOM = ./src/DIATOM/
OBJ_PATH_DIATOM = ./obj/DIATOM/
INC_PATH_DIATOM = -I ./include/DIATOM/

$(OBJ_PATH_DIATOM)matrix.o: $(SRC_PATH_DIATOM)matrix.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ $(INC_PATH_DIATOM) $(INC_PATH_DIP) $(INC_PATH_POT) 
######################################################################

######################################################################
# Dipole moment
#
SRC_PATH_DIP = ./src/DIP/
OBJ_PATH_DIP = ./obj/DIP/
INC_PATH_DIP = -I ./include/DIP/

$(OBJ_PATH_DIP)co2_ar_dipole.o: $(SRC_PATH_DIP)co2_ar_dipole.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_DIP)

$(OBJ_PATH_DIP)ar_he_dip.o: $(SRC_PATH_DIP)ar_he_dip.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_DIP)

$(OBJ_PATH_DIP)ar_he_dip_buryak_fit.o: $(SRC_PATH_DIP)ar_he_dip_buryak_fit.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_DIP)

$(OBJ_PATH_DIP)ar_he_dip_derivative.o: $(SRC_PATH_DIP)ar_he_dip_derivative.cpp
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

$(OBJ_PATH_MISC)file.o: $(SRC_PATH_MISC)file.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_MISC)

$(OBJ_PATH_MISC)parameters.o: $(SRC_PATH_MISC)parameters.cpp ./include/MISC/parameters.h 
	$(CCXX) $(CXXFLAGS) -c $< -o $@ $(INC_PATH_MISC)

$(OBJ_PATH_MISC)spectrum_info.o: $(SRC_PATH_MISC)spectrum_info.cpp
	$(MPICC) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH_MISC)
######################################################################

######################################################################
# path to miscellaneous functions
MISC_OBJ = aux.o fft.o 
MISC_OBJECTS = $(patsubst %,$(OBJ_PATH_MISC)%,$(MISC_OBJ))

MISC_OBJ_CUT = fft.o
MISC_OBJECTS_CUT = $(patsubst %,$(OBJ_PATH_MISC)%,$(MISC_OBJ_CUT))
######################################################################

######################################################################
# CO2-Ar trajectory spectrum
#
POT_OBJ = psp_pes.o leg_arr.o 
POT_OBJECTS = $(patsubst %,$(OBJ_PATH_POT)%,$(POT_OBJ))

DIP_OBJ = co2_ar_dipole.o 
DIP_OBJECTS = $(patsubst %,$(OBJ_PATH_DIP)%,$(DIP_OBJ))

CO2AR_OBJ = matrix.o 
CO2AR_OBJECTS = $(patsubst %,$(OBJ_PATH_CO2AR)%,$(CO2AR_OBJ))

CO2AR_TRAJ_OBJ = $(GEAR_OBJECTS) $(POT_OBJECTS) $(DIP_OBJECTS) $(CO2AR_OBJECTS) $(MISC_OBJECTS) $(CONSTANTS_OBJECTS)

CO2AR_TRAJ_INC = $(INC_PATH_GEAR) $(INC_PATH_MISC) $(INC_PATH_CO2AR) $(INC_PATH_CONSTANTS) $(INC_PATH_POT) $(INC_PATH_DIP)

co2ar_trajs: $(CO2AR_TRAJ_OBJ) $(SRC_PATH_CO2AR)main.cpp
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $^ -o $@ $(CO2AR_TRAJ_INC) 

# with euler angles
CO2AR_OBJ_EULER = matrix_euler.o
CO2AR_OBJECTS_EULER = $(patsubst %,$(OBJ_PATH_CO2AR)%,$(CO2AR_OBJ_EULER))

CO2AR_TRAJ_OBJ_EULER = $(GEAR_OBJECTS_DIATOM) $(POT_OBJECTS) $(DIP_OBJECTS) $(CO2AR_OBJECTS_EULER) $(MISC_OBJECTS) $(CONSTANTS_OBJECTS)

co2ar_trajs_euler: $(CO2AR_TRAJ_OBJ_EULER) $(SRC_PATH_CO2AR)main_euler.cpp
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $^ -o $@ $(CO2AR_TRAJ_INC)
######################################################################

######################################################################
# Diatomic trajectory spectrum
#
DIATOM_POT_OBJ = ar_he_pes.o ar_he_pes_derivative.o
DIATOM_POT_OBJECTS = $(patsubst %,$(OBJ_PATH_POT)%,$(DIATOM_POT_OBJ))

DIATOM_DIP_OBJ = ar_he_dip_buryak_fit.o 
DIATOM_DIP_OBJECTS = $(patsubst %,$(OBJ_PATH_DIP)%,$(DIATOM_DIP_OBJ))

DIATOM_OBJ = matrix.o
DIATOM_OBJECTS = $(patsubst %,$(OBJ_PATH_DIATOM)%,$(DIATOM_OBJ))

# using special Gear method without additional parameters 
DIATOM_TRAJ_OBJ = $(GEAR_OBJECTS_DIATOM) $(DIATOM_POT_OBJECTS) $(DIATOM_DIP_OBJECTS) $(DIATOM_OBJECTS) $(MISC_OBJECTS) 
DIATOM_TRAJ_INC = $(INC_PATH_POT) $(INC_PATH_DIP) $(INC_PATH_GEAR) $(INC_PATH_MISC) $(INC_PATH_DIATOM) $(INC_PATH_CONSTANTS)

diatom_trajs: $(DIATOM_TRAJ_OBJ) $(SRC_PATH_DIATOM)main.cpp
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $^ -o $@ $(DIATOM_TRAJ_INC)

buryak: $(DIATOM_TRAJ_OBJ) $(SRC_PATH_DIATOM)main_buryak.cpp
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $^ -o $@ $(DIATOM_TRAJ_INC)

formatted_coords: $(DIATOM_TRAJ_OBJ) $(SRC_PATH_DIATOM)formatted_coords.cpp
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $(INCLUDE_GSL) $^ -o $@ $(DIATOM_TRAJ_INC)

SIMPLE_MC_OBJ = $(OBJ_PATH_MISC)file.o $(OBJ_PATH_MISC)parameters.o $(OBJ_PATH_MISC)spectrum_info.o 
simple_monte_carlo: $(DIATOM_TRAJ_OBJ) $(SIMPLE_MC_OBJ) $(SRC_PATH_DIATOM)simple_monte_carlo.cpp 
	$(MPICC) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_FFTW) $^ -o $@ $(DIATOM_TRAJ_INC)
######################################################################

######################################################################
# plot spectral function and binary adsorption coeff
#
make_spectrum: $(SRC_PATH_MISC)make_spectrum.cpp
	$(CCXX) $(CXXFLAGS_NOOPT) $(INCLUDE_GNUPLOT_IO) $(LINK_BOOST) $^ -o $@ $(INC_PATH_CONSTANTS)
######################################################################

######################################################################
visualize: $(SRC_PATH_DIATOM)visualize_trajectory.cpp
	$(CCXX) $(CXXFLAGS_NOOPT) $(OPENGL_LFLAGS) $^ -o $@
######################################################################

######################################################################
SRC_PATH_TEST = ./tests/

potential_test: $(SRC_PATH_TEST)potential_test.cpp $(DIATOM_DIP_OBJECTS) $(DIATOM_POT_OBJECTS) 
	$(CCXX) $(CXXFLAGS_NOOPT) $^ -o $@ $(INC_PATH_DIP) $(INC_PATH_POT)

fft_test: $(SRC_PATH_TEST)fft_test.cpp $(MISC_OBJECTS)
	$(CCXX) $(CXXFLAGS_NOOPT) $^ -o $@ $(INC_PATH_MISC) $(INCLUDE_FFTW) $(INCLUDE_GSL) $(INCLUDE_GNUPLOT_IO) $(LINK_BOOST)

real_fft_test: $(SRC_PATH_TEST)real_fft_test.cpp $(MISC_OBJECTS)
	$(CCXX) $(CXXFLAGS_NOOPT) $^ -o $@ $(INC_PATH_MISC) $(INCLUDE_FFTW) $(INCLUDE_GSL)
######################################################################




