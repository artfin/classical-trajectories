#pragma once

#include <iostream>
#include "matrix.h"

#include <vector>
#include "parameters.h"

// Gear header files
#include "basis.h"
#include "vmblock.h"
#include "gear.h"
#include "t_dgls.h"

class Trajectory
{
public:
	int N;
	Parameters parameters;

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
	
	// dipole moment in laboratory frame
	std::vector<double> temp{ std::vector<double>(3) };
	std::vector<double> dipx;
	std::vector<double> dipy;
	std::vector<double> dipz;

	void syst( REAL t, REAL *y, REAL *f );
	void run_trajectory( void );
	Trajectory ( Parameters& parameters, const int& N );
};
