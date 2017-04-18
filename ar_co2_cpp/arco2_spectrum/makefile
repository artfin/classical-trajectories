all: trajectories.exe clean

CXXFLAGS := -O3 -march=native -mtune=native

trajectories.exe: trajectories.o basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o f.o  leg_arr.o matrix.o
	g++ $(CXXFLAGS) trajectories.o basis_r.o awp.o vmblock.o fgauss.o t_dgls.o psp_pes.o gear.o f.o  leg_arr.o matrix.o -o trajectories.exe


awp.o: 
	g++ $(CXXFLAGS) -c ../jean-pierre/awp.cpp

basis_r.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/basis_r.cpp 

gear.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/gear.cpp

t_dgls.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/t_dgls.cpp

vmblock.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/vmblock.cpp

fgauss.o: 
	g++ $(CXXFLAGS) -c  ../jean-pierre/fgauss.cpp

leg_arr.o: 
	g++ $(CXXFLAGS) -c  ../parker_snow/leg_arr.cpp

psp_pes.o: 
	g++ $(CXXFLAGS) -c  ../parker_snow/psp_pes.cpp

f.o: f.cpp
	g++ $(CXXFLAGS) -c f.cpp

matrix.o: matrix.cpp
	g++ $(CXXFLAGS) -c matrix.cpp


trajectories.o: trajectories.cpp
	g++ $(CXXFLAGS) -c trajectories.cpp

clean:
	rm *.o
