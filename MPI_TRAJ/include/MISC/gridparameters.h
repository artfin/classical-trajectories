#pragma once

#include <iostream>

using std::cout;
using std::endl;

class GridParameters
{
	public:
		double B_MIN;
		double B_MAX;

		double V0_MIN;
		double V0_MAX;

		GridParameters( );
		~GridParameters( );
};
