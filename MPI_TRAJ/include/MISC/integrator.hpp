#pragma once

#include "hep/mc-mpi.hpp"
#include "integrand.hpp"

#include <iostream>

using std::cout;
using std::endl;

class Integrator
{
public:
	Integrand& integrand;
	int world_rank;

	Integrator( Integrand& integrand, const int& world_rank ) : integrand(integrand), world_rank(world_rank) 
	{
	}

	void set_callback( void )
	{
		hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);
	}

	double run_integration( const int& niter, const int& ndots )
	{
		auto results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>( integrand, integrand.DIM ),
			std::vector<std::size_t>( niter, ndots )
		);

		auto result = hep::cumulative_result0( results.begin() + 1, results.end() );

		//if ( world_rank == 0 )
		//{
			//cout << ">> value: " << result.value() << "; error: " << result.error() << "; relative error (%): " << (double) result.error() / result.value() * 100 << endl;
		//}

		return result.value();
	}
};
