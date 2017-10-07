CCXX = g++
CXXFLAGS = -std=c++11 -lfftw3 -lm

dft: dft.cpp
	$(CCXX) $(CXXFLAGS) dft.cpp -o $@

fft-test: fft-test.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@
