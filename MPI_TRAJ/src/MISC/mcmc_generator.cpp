#include "mcmc_generator.hpp"

MCMC_generator::MCMC_generator( std::function<double(VectorXd)> f, 
				 				const int& DIM, 
								const double& alpha, 
								const int& subchain_length, 
								vector<pair<int, double>>& to_wrap ) : 
		DIM(DIM), 
		alpha(alpha), 
		subchain_length(subchain_length),
		to_wrap(to_wrap),
		f(f)
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

void MCMC_generator::nextGaussianVec( VectorXd& v, VectorXd& mean )
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

VectorXd MCMC_generator::metro_step( VectorXd& x )
{
	VectorXd prop( DIM );
	nextGaussianVec( prop, x );

	if ( nextDouble(0.0, 1.0) < std::min( 1.0, f(prop) / f(x) ))
	{
		return prop;
	}

	return x;
}

void MCMC_generator::burnin( VectorXd initial_point, const int& burnin_length )
{
	VectorXd x = metro_step( initial_point );

	for ( size_t i = 0; i < burnin_length; i++ )
	{
		x = metro_step( x );
	}

	current_point = x;
	burnin_done = true;
}

VectorXd MCMC_generator::generate_point( )
{
	if ( burnin_done == false )
	{
		std::cout << "Burnin is not done!" << std::endl;
		exit( 1 );
	}

	int moves = 0;
	VectorXd x = current_point;
	VectorXd xnew;

	while ( moves < subchain_length )
	{
		xnew = metro_step( x );

		if ( to_wrap != DEFAULT_VECTOR )
		{
			int curr_var; // number of current variable
			double curr_max; // maximum of current variable

			for ( size_t i = 0; i < to_wrap.size(); i++ )
			{
				curr_var = to_wrap[i].first;
				curr_max = to_wrap[i].second;

				xnew(curr_var) = wrapMax( xnew(curr_var), curr_max );
			}
		}

		if ( xnew != x )
		{
			x = xnew;
			moves++;
		}
	}	

	current_point = x;
	return x;
}

