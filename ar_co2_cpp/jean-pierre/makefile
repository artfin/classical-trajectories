all: gear_solver.exe clean

#CXXFLAGS := -Wall

#gear_solver.exe: awp.o basis_r.o gear.o t_dgls.o vmblock.o fgauss.o mgear.o
#	g++ $(CXXFLAGS) awp.o basis_r.o gear.o t_dgls.o vmblock.o fgauss.o mgear.o -o gear_solver.exe


gear_solver.exe: awp.o basis_r.o gear.o t_dgls.o vmblock.o fgauss.o  system.o f.o psp_pes_tapenade_d.o
	g++ $(CXXFLAGS) awp.o basis_r.o gear.o t_dgls.o vmblock.o fgauss.o system.o psp_pes_tapenade_d.o f.o -o gear_solver.exe


awp.o: awp.cpp
	g++ $(CXXFLAGS) -c awp.cpp

basis_r.o: basis_r.cpp
	g++ $(CXXFLAGS) -c basis_r.cpp 

gear.o: gear.cpp
	g++ $(CXXFLAGS) -c gear.cpp

t_dgls.o: t_dgls.cpp
	g++ $(CXXFLAGS) -c t_dgls.cpp

vmblock.o: vmblock.cpp
	g++ $(CXXFLAGS) -c vmblock.cpp

fgauss.o: fgauss.cpp
	g++ $(CXXFLAGS) -c fgauss.cpp

mgear.o: mgear.cpp	
	g++ $(CXXFLAGS) -c mgear.cpp

system.o: system.cpp	
	g++ $(CXXFLAGS) -c system.cpp

f.o: 
	g++ $(CXXFLAGS) -c sys\f.cpp

psp_pes_tapenade_d.o: 
	g++ $(CXXFLAGS) -c sys\psp_pes_tapenade_d.cpp



clean:
	rm *.o
