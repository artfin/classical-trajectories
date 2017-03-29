CCXX = g++
CXXFLAGS = -std=c++11 -I/usr/local/include/ -L/usr/local/lib/
DEPS = mpfr.h

LIBS = -lsymengine -lgmp -lgmpxx -lsymengine -lz 

%.o: %.cpp $(DEPS)
	$(CCXX) -c -o $@ $< $(CXXFLAGS)

main: main.o
	g++ -o main main.o $(CXXFLAGS) $(LIBS)
	# $(CCXX) $^ -o $@ $(CXXFLAGS) $(LIBS)
