#include "ab_initio_potential.h"
#include "hep/mc.hpp"

#include <math.h> 
#include <cstddef>
#include <iostream>
#include <vector>

// computes power of any integer to any POSITIVE power
static long long int integer_pow(int x, int n)
{
	long long int r = 1;
	while ( n-- )
	{
		r *= x;
	}

	return r;
}

// hartree to joules
const double HTOJ = 4.35974417 * pow(10, -18);
// boltzmann constant
const long double BOLTZCONST = 1.38064852 * pow(10, -23);
// current temperature
const double Temperature = 200;

double integrand(hep::mc_point<double> const& x)
{
	double potential_value = ab_initio_pot(x.point()[0], x.point()[1]);

	return ( 1 - exp( - potential_value / (BOLTZCONST * Temperature) )) * pow(x.point()[0], 2) * sin(x.point()[1]);
}


int main()
{
	std::cout << ">> computing integral of AB INITIO potential" << std::endl;
	
	std::cout << ">> Just checking that some constants are normally introduced.." << std::endl;
	std::cout << ">> HTOJ = " << HTOJ << std::endl;
	std::cout << ">> BOLTZCONST = " << BOLTZCONST << std::endl;
	std::cout << "----------------------------------" << std::endl;

	// set the verbose vegas callback function
	hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);

	// perform 5 iteration with 1000 calls each
	// this function also calls vegas_verbose_callback after each iteration which in turn prints the individual iterations
	auto results = hep::vegas(
		hep::make_integrand<double>(integrand, 2),
		std::vector<std::size_t>(5, 1000)
	);

	// results contains the estimations for each iteration. We could take the
	// result from last iteration, but here we instead choose to combine the
	// results of all iterations but the first one in a cumulative
	auto result = hep::cumulative_result0(results.begin() + 1, results.end());
	double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

	// print the cumulative result
	std::cout << ">> cumulative result (excluding the first iteration): \n>> N=" << result.calls() << " I=" << result.value() << " +- " << result.error() << " chi^2/dof=" << chi_square_dof << std::endl;

	return 0;
}


