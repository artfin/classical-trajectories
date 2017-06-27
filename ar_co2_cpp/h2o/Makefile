all: one_dim_water

CCXX = g++
CXXFLAGS = -I /usr/local/include/eigen3 -std=c++11

one_dim_water: matrix.o
	$(CCXX) $(CXXFLAGS) matrix.o -o $@

matrix.o: matrix.cpp
	$(CCXX) $(CXXFLAGS) -c matrix.cpp

clean: 
	rm *.o
