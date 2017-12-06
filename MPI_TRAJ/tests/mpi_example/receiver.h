#pragma once

#include <mpi.h>
#include <iostream>

class Receiver
{
public:
	int x;
	Receiver( ) { };
	void receive( void );
};

