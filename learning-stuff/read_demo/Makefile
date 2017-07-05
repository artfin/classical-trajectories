EXECS = reader boost1 boost2 boost3 boost4 boost5 boost6 alex
CCXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O3

all: ${EXECS}

reader: reader.cc
	${CCXX} ${CXXFLAGS} -o reader reader.cc

boost1: boost1.cc
	${CCXX} ${CXXFLAGS} -o boost1 boost1.cc

boost2: boost2.cc
	${CCXX} ${CXXFLAGS} -o boost2 boost2.cc

boost3: boost3.cc
	${CCXX} ${CXXFLAGS} -o boost3 boost3.cc

boost4: boost4.cc
	${CCXX} ${CXXFLAGS} -o boost4 boost4.cc

boost5: boost5.cc
	${CCXX} ${CXXFLAGS} -o boost5 boost5.cc

boost6: boost6.cc
	${CCXX} ${CXXFLAGS} -o boost6 boost6.cc

alex: alex.cc
	${CCXX} ${CXXFLAGS} -o alex alex.cc

clean:
	rm -f ${EXECS}

