#pragma once

#include <mpi.h>
#include <iostream>

class Sender
{
public:
	int x;
	
	Sender( int x ) : x(x) { };

	void send( void );
};


