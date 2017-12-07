#pragma once

#include <mpi.h>

#include <iostream>
#include <vector>
#include <string>
#include <functional>

class SpectrumInfo
{
public:
	SpectrumInfo( int size ) : size(size)
	{
		reserve_space();
	}	
	
	SpectrumInfo( std::function<double(double, double)> corrector ) : corrector(corrector) { }	

	SpectrumInfo( ) { }

	int size;
	std::function<double(double, double)> corrector;

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

	void zero_out_chunk( void );
	void zero_out_package( void );
	void add_chunk_to_total_and_zero_out( void );

	void correct( SpectrumInfo& classical, double omega, double freq_step, double kT, bool beware_zero = false );

	void send( void );
	double receive( int source = 0, bool define_source = false );
};


