#pragma once

#include <iostream>
#include <random>
#include <string>

using std::cout;
using std::endl;
using std::string;

static std::mt19937 generator;

struct ICPoint
{
	int counter;
	double b;
	double v0;
};

class Parameters
{
	public:
		double B_MIN;
		double B_MAX;
		int B_PARTS;

		double V0_MIN;
		double V0_MAX;
		int V0_PARTS;

		int CYCLES;
		int CYCLE_POINTS;

		double STDDEV_MAX;

		string specfunc_filename = "";
		string spectrum_filename = "";
		string m2_filename = "";

		Parameters( );
		~Parameters( );

		void show_parameters( void );
		
		double nextDouble( const double &min, const double &max );	
		ICPoint generate_uniform_point( void );
		ICPoint generate_uniform_point( const double& bmin,
										const double& bmax,
										const double& v0min,
										const double& v0max );
};
