#pragma once

#include <mpi.h>

#include <iostream>
#include <vector>
#include <string>

class SpectrumInfo
{
public:
	SpectrumInfo( int size ) : size(size)
	{
		reserve_space();
	}		

	int size;

	std::vector<double> specfunc_package;
	std::vector<double> spectrum_package;

	std::vector<double> specfunc_chunk;
	std::vector<double> spectrum_chunk;

	std::vector<double> specfunc_total;
	std::vector<double> spectrum_total;

	double m2_package;
	double m2_chunk;
	double m2_total;

	void reserve_space( void );

	void multiply_chunk( double multiplier );
	void add_chunk_to_total( void );
	void add_package_to_chunk( void );

	void zero_out( std::vector<double>& v );
	void zero_out_chunk( void );
	void add_chunk_to_total_and_zero_out( void );

	double receive( int source = -1 );
};


