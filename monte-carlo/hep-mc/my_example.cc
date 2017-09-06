#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

// the function that shall be integrated
double func(hep::mc_point<double> const& x)
{
	double x_new = std::tan(M_PI / 2 * x.point()[0]);

	return M_PI / 2 * (1 + std::pow(x_new, 2)) * std::exp( -x_new ); 
}

int main()
{
	std::cout << ">> computing integral of exp(-x) from 0 to +inf\n";
	
	hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);

	auto results = hep::vegas(
		hep::make_integrand<double>(func, 1),
		std::vector<std::size_t>(5, 10000)
	);

	auto result = hep::cumulative_result0(results.begin() + 1, results.end());
	double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

	std::cout << ">> cumulative result (excluding the first iteration):" << std::endl;
	std::cout << ">> N=" << result.calls() << std::endl;
	std::cout << ">> I=" << result.value() << "+-" << result.error() << std::endl;
	std::cout <<">> chi^2/dof=" << chi_square_dof << std::endl;

	return 0;
}
