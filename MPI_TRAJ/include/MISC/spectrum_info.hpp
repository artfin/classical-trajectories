#pragma once

#include <mpi.h>

#include <iostream>
#include <vector>
#include <string>
#include <functional>

// checking if dir exists
#include <dirent.h>
#include <sys/stat.h> // mkdir

#include <assert.h>
#include <fstream>

#include "parameters.h"

class SpectrumInfo
{
public:
	SpectrumInfo( int size, string modifier ) : size(size), modifier(modifier)
	{
		reserve_space();
	}	
	
	SpectrumInfo( std::function<double(double, double)> corrector ) : corrector(corrector) { }	

	SpectrumInfo( ) { } 

	int size;
	std::string modifier; // classical, d1, d2, d3 
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
	void multiply_total( const double& multiplier );
	void add_chunk_to_total( void );
	void add_package_to_chunk( void );
	void add_package_total( void );

	void zero_out_chunk( void );
	void clear_package( void );

	void correct( SpectrumInfo& classical, double omega, double freq_step, double kT, bool beware_zero = false );

	void send( void );
	void receive( int& source, bool define_source = false );

	void save( std::vector<double>& v1, std::vector<double>& v2, std::string filename );
	void save( const double& m2, std::string filename );

	void saving_procedure( Parameters& parameters, std::vector<double>& freqs, std::string filename = "" );
	std::string modify_filename( std::string filename, std::string modifier );
	bool check_dir_exists( std::string dirname );
};


