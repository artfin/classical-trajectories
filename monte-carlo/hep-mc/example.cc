#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

// the function that shall be integrated
double square(hep::mc_point<double> const& x)
{
	return 3.0 * x.point()[0] * x.point()[0];
}

int main()
{
	std::cout << ">> computing integral of 3*x^2 from 0 to 1 which is 1\n";

	return 0;	
}
