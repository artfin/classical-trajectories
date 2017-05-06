CCXX = g++
CXXFLAGS = -I /usr/local/include/eigen3/ -std=c++11

ham: hamiltonian.cpp
	$(CCXX) -O3 -march=native -mtune=native -funroll-loops -o $@ $^ $(CXXFLAGS)

