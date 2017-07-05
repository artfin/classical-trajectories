.PHONY: all clean

EXECS = main

CCXX = g++
CXXFLAGS = -std=c++11 -I/home/artfin/Downloads/hep-mc-0.5/include -lm

all: ${EXECS}

main: ab_initio_potential.o main.cc
	$(CCXX) $(CXXFLAGS) $^ -o $@

potential.o: ab_initio_potential.c
	$(CCXX) $(CXXFLAGS) -c $@ -o $@

clean: 
	rm -rf *.o
