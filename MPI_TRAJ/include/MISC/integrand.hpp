#pragma once

#include "hep/mc.hpp"
#include "limits.hpp"

#include <iostream>
#include <vector>
#include <functional>

#include <Eigen/Dense>

using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::function;
using std::vector;

class Integrand
{
public:
	int DIM;
	function<double(VectorXd)> f;
	Limits limits;

	Integrand( function<double(VectorXd)> f, const int& DIM ) :
		f(f), DIM(DIM) { }

	Limits* set_limits( void ) { return &limits; }

	double operator()( hep::mc_point<double> const& x )
	{
		// x transformed from [0, 1]^n cube to normal coordinates
		VectorXd xtr( DIM );
		
		double jac = 1;

		for ( size_t i = 0; i < DIM; i++ )
		{
			if ( limits.limits[i].chvarType == Limit::chvarTypes::INFINF )
			{
				xtr(i) = tan( M_PI * (x.point()[i] - 0.5) );
			   	jac *= M_PI * (1 + xtr(i) * xtr(i));	
			}
			
			if ( limits.limits[i].chvarType == Limit::chvarTypes::FINITE )
			{
				xtr(i) = ( limits.limits[i].ub - limits.limits[i].lb ) * x.point()[i];
				jac *= ( limits.limits[i].ub - limits.limits[i].lb );
			}

			if ( limits.limits[i].chvarType == Limit::chvarTypes::ZEROINF )
			{
				xtr(i) = tan( M_PI / 2 * x.point()[i] );
				jac *= M_PI / 2.0 * ( 1 + xtr(i) * xtr(i) );
			}
		}

		return f( xtr ) * jac; 
	}
};
