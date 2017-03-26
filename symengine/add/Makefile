# Makefile 1 
# main: main.cpp add.cpp
#
# It is not very efficient because it will recompile all files ignoring that some of them haven't changed
#
#
# Makefile 2

# macro CC is the C compiler to use, and CFLAGS is the list of flags to pass to the compilation command
# by putting the object file -- main.o and add.o in the dependency list and in the rule, 
# make knows it must first compile the .cpp versions individually, and then build the executable main
#CC=g++
#CFLAGS=-I.

#main: main.o add.o
	#$(CC) -o main main.o add.o -I.

# But we need to all .cpp files to be dependent on certain .h files.

# Makefile 3

CC = g++
CFLAGS = -std=c++11
DEPS = add.h

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: main.o add.o
	g++ -o main main.o add.o -I.


