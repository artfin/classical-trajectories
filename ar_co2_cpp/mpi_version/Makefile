EXECS = gather send_recv 
MPICC?=mpicc
CXXFLAGS = -std=c++0x -lstdc++

all: ${EXECS}

gather: gather.cc
	${MPICC} ${CXXFLAGS} -o gather gather.cc 

send_recv: send_recv.cc
	${MPICC} ${CXXFLAGS} -o send_recv send_recv.cc

clean:
	rm -f ${EXECS}
