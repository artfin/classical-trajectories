#pragma once

#include <iostream>
#include <random>
#include <string>

using std::cout;
using std::endl;
using std::string;

static std::mt19937 generator( 27717 );

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

		int NPOINTS;

		double RDIST;
		double sampling_time;
		int MaxTrajectoryLength;
		double FREQ_MAX;
		bool use_S_matrix;

		bool d1_status;
		bool d2_status;
		bool d3_status;
		bool d4_status;

		double Temperature;
			
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
