.PHONY: all clean

EXECS = main

# MPI compiler wrapper
MPICC = mpic++ 

# usual compiler
CCXX = g++

CXXFLAGS = -std=c++11 -lstdc++ -lfftw3 -I /usr/local/include/eigen3 -march=native -mtune=native

all: ${EXECS}

objects = basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o f.o leg_arr.o matrix.o co2_ar_dipole.o

main: $(objects) 
	$(MPICC) $(objects) $(CXXFLAGS) -c main.cc -o main 

co2_ar_dipole.o:
	$(CCXX) $(CXXFLAGS) -c ../arco2_spectrum/co2_ar_dipole.cpp

matrix.o:
	$(CCXX) $(CXXFLAGS) -c ../arco2_spectrum/matrix.cpp

awp.o:
	$(CCXX) $(CXXFLAGS) -c ../jean-pierre/awp.cpp 

basis_r.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/basis_r.cpp 

gear.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/gear.cpp

t_dgls.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/t_dgls.cpp

vmblock.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/vmblock.cpp

fgauss.o: 
	$(CCXX) $(CXXFLAGS) -c  ../jean-pierre/fgauss.cpp

leg_arr.o: 
	$(CCXX) $(CXXFLAGS) -c ../parker_snow/leg_arr.cpp 

psp_pes.o: 
	$(CCXX) $(CXXFLAGS) -c  ../parker_snow/psp_pes.cpp

f.o: 
	$(CCXX) $(CXXFLAGS) -c ../arco2_spectrum/f.cpp 

clean:
	rm -f ${EXECS}

