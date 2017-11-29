#include "parameters.h"

Parameters::Parameters( )
{
	//cout << "Initialized GridParameters class" << endl;
}

Parameters::~Parameters()
{
	//cout << "Destructing GridParameters class" << endl;
}

void Parameters::show_parameters( void )
{
	cout << "Monte-Carlo parameters:" << endl;
	cout << "CYCLES: " << this->CYCLES << endl;
	cout << "CYCLE_POINTS: " << this->CYCLE_POINTS << endl;
	cout << "STDDEV_MAX: " << this->STDDEV_MAX << endl;
   	cout << "#################################" << endl;	
	cout << "Gridparameters: " << endl;
	cout << "B_MIN: " << this->B_MIN << endl;
	cout << "B_MAX: " << this->B_MAX << endl;
	cout << "V0_MIN: " << this->V0_MIN << endl;
	cout << "V0_MAX: " << this->V0_MAX << endl;
	cout << "################################" << endl;
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
