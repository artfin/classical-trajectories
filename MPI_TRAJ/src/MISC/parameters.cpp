#include "parameters.h"

Parameters::Parameters( )
{
	// constructor
}

Parameters::~Parameters()
{
	// destructor
}

void Parameters::show_parameters( void )
{
   	cout << "#######################" << endl;	
   	cout << "Conditions parameters: " << endl;
	cout << "Temperature: " << this->Temperature << endl;
	cout << "#######################" << endl;	
	cout << "Trajectory parameters: " << endl;
	cout << "RDIST (initial distance): " << this->RDIST << endl;
	cout << "sampling time: " << this->sampling_time << endl;
	cout << "MaxTrajectoryLength: " << this->MaxTrajectoryLength << endl;
	cout << "FREQ_MAX (maximum frequency in spectrum): " << this->FREQ_MAX << endl;
	cout << "#######################" << endl;	
	cout << "Monte-Carlo parameters:" << endl;
	cout << "NPOINTS" << this->NPOINTS << endl;
   	cout << "#######################" << endl;	
	cout << "Gridparameters: " << endl;
	cout << "B_MIN: " << this->B_MIN << endl;
	cout << "B_MAX: " << this->B_MAX << endl;
	cout << "B_PARTS: " << this->B_PARTS << endl;
	cout << "V0_MIN: " << this->V0_MIN << endl;
	cout << "V0_MAX: " << this->V0_MAX << endl;
	cout << "V0_PARTS: " << this->V0_PARTS << endl;
	cout << "########################" << endl;
	cout << "Files parameters: " << endl;
	cout << "specfunc filename: " << this->specfunc_filename << endl;
	cout << "spectrum filename: " << this->spectrum_filename << endl;
	cout << "m2 filename: " << this->m2_filename << endl;
}


double Parameters::nextDouble( const double &min, const double &max )
{
	std::uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}


ICPoint Parameters::generate_uniform_point( void )
{
	ICPoint p;

	p.b = nextDouble( this->B_MIN, this->B_MAX );
	p.v0 = nextDouble( this->V0_MIN, this->V0_MAX );

	return p;
}

ICPoint Parameters::generate_uniform_point( const double& bmin,
			   								const double& bmax,
											const double& v0min,
											const double& v0max )
{
	ICPoint p;

	p.b = nextDouble( bmin, bmax );
	p.v0 = nextDouble( v0min, v0max );

	return p;
}


