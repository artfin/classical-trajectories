CCXX = g++
CXXFLAGS = -I /usr/local/include/eigen3/ -std=c++11

test: hamiltonian.cpp
	$(CCXX) -O3 -march=native -mtune=native -funroll-loops -o $@ $^ $(CXXFLAGS)
