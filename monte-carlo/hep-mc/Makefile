EXECS = introduction introduction_mpi stop_after_precision vegas_write_pdf 

MPICC = mpic++
CCXX = g++
CXXFLAGS = -std=c++11 -I/home/artfin/Downloads/hep-mc-0.5/include 

all: ${EXECS}

introduction: introduction.cc
	$(CCXX) $(CXXFLAGS) $^ -o $@

introduction_mpi: introduction_mpi.cc
	$(MPICC) $(CXXFLAGS) $^ -o $@

stop_after_precision: stop_after_precision.cc
	$(CCXX) $(CXXFLAGS) $^ -o $@

vegas_write_pdf: vegas_write_pdf.cc
	$(CCXX) $(CXXFLAGS) $^ -o $@
