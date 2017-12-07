#include "trajectory.hpp"

Trajectory::Trajectory( int N )
{
	this->N = N;

	
	vmblock = vminit();
	y0 = (REAL*) vmalloc(vmblock, VEKTOR, N, 0);
}

void Trajectory::run_trajectory()
{
	// out of memory?
	if ( !vmcomplete(vmblock) )
	{ 
		cout << "mgear: out of memory" << endl;
		return;
	}
	
}
