CCXX = g++
CXXFLAGS = -std=c++11 -lm
STFLAGS = -I /usr/local/include/eigen3

all: ziggurat rand_rot

ziggurat: ziggurat.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

rand_rot: rand_rot.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@
