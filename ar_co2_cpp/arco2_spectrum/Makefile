all: trajectories clean

CCXX = g++
CXXFLAGS = -I /usr/local/include/eigen3/ -std=c++11 -O3 -march=native -mtune=native

trajectories: basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o f.o  leg_arr.o matrix.o co2_ar_dipole.o trajectories.o
	$(CCXX) $(CXXFLAGS) trajectories.o basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o f.o  leg_arr.o matrix.o  co2_ar_dipole.o -o $@ 


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
	$(CCXX) $(CXXFLAGS) -c  ../parker_snow/leg_arr.cpp

psp_pes.o: 
	$(CCXX) $(CXXFLAGS) -c  ../parker_snow/psp_pes.cpp

f.o: f.cpp
	$(CCXX) $(CXXFLAGS) -c f.cpp

co2_ar_dipole.o: 
	$(CCXX) $(CXXFLAGS) -c co2_ar_dipole.cpp

matrix.o: matrix.cpp
	$(CCXX) $(CXXFLAGS) -c matrix.cpp

trajectories.o: trajectories.cpp
	$(CCXX) $(CXXFLAGS) -c trajectories.cpp

clean:
	rm *.o
