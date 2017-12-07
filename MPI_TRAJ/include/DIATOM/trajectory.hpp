#pragma once

#include <iostream>

// Gear header files
#include "basis.h"
#include "vmblock.h"
#include "gear.h"
#include "t_dgls.h"

class Trajectory
{
public:
	int N;
	
	REAL epsabs = 1E-13;    //  absolute error bound
	REAL epsrel = 1E-13;    //  relative error bound    
	REAL t0;        // left edge of integration interval
	REAL *y0;       // [0..n-1]-vector: initial value, approxim. 
	REAL h;         // initial, final step size
	REAL xend;      // right edge of integration interval 

	long fmax = 1e8;      // maximal number of calls of right side in gear4()
	long aufrufe;   // actual number of function calls
	int  fehler;    // error code from umleiten(), gear4()

	void *vmblock;  // List of dynamically allocated vectors
	

	Trajectory ( int N );
};
