CCXX = g++
CXXFLAGS = -std=c++11 -lm -O3
STFLAGS = -I /usr/local/include/eigen3

all: ziggurat rand_rot diatomics_distr

ziggurat: ziggurat.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

rand_rot: rand_rot.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $^ -o $@

diatomics_distr: diatomics_distr.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $^ -o $@

