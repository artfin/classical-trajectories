#include "mcmc_generator.hpp"

MCMC_generator::MCMC_generator( function<double(VectorXf)> f, const int& DIM, const double& alpha, const int& n ) : DIM(DIM), alpha(alpha), f(f), n(n)
{
}

MCMC_generator::~MCMC_generator()
{
}

double MCMC_generator::nextDouble( const double& min, const double& max )
{
	std::uniform_real_distribution<double> distributuion( min, max );
	return distributuion( this->generator );
}

void MCMC_generator::nextGaussianVec( VectorXf &v, VectorXf& mean )
{
	for ( size_t i = 0; i < DIM; i++ )
	{
		std::normal_distribution<double> d( mean(i), alpha );
		v(i) = d( generator );
	}
}

double MCMC_generator::wrapMax( const double& x, const double& max )
{
	return fmod( max + fmod(x, max), max );
}

VectorXf MCMC_generator::metro_step( VectorXf& x )
{
	VectorXf prop( DIM );
	nextGaussianVec( prop, x );

	if ( nextDouble(0.0, 1.0) < std::min( 1.0, f(prop) / f(x) ))
	{
		return prop;
	}

	return x;
}

VectorXf MCMC_generator::generate_point( vector<pair<int, double>>& to_wrap )
{
	
}
