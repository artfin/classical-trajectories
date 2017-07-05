#include "hep/mc-mpi.hpp"

#include <mpi.h>

#include <cstddef>
#include <iostream>
#include <vector>

// the function that shall be integrated
double square(hep::mc_point<double> const& x)
{
	return 3.0 * x.point()[0] * x.point()[0];
}

int main(int argc, char* argv[])
{
	// Initialize MPI
	MPI_Init(&argc, &argv);

	// obtaining the rank of the current process
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// printing info
	if ( rank == 0 )
	{
		std::cout << ">> computing integral of 3*x^2 from 0 to 1 which is 1.0" << std::endl;
	}

	// set the verbose vegas callback function
	hep::mpi_vegas_callback<double>(hep::mpi_vegas_verbose_callback<double>);

	// perform 5 iteration with 1000 calls each
	auto results = hep::mpi_vegas(
		MPI_COMM_WORLD,
		hep::make_integrand<double>(square, 1),
		std::vector<std::size_t>(5, 10000000)
	);

	// results contains estimations for each iteration. calculating cumulation w/out 1 result
	auto result = hep::cumulative_result0(results.begin() + 1, results.end());
	double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

	// print the cumulative result
	if ( rank == 0 )
	{
		std::cout << ">> cumulative result (excluding first iteration):\n>> N=" << result.calls() << " I=" << result.value() << " +- " << result.error() << " chi^2/dof=" << chi_square_dof << std::endl;
	}

	// MPI clean up
	MPI_Finalize();
	
	return 0;
}
