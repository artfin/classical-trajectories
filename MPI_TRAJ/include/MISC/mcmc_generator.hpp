#pragma once

#include <iostream>
#include <random>
#include <functional>
#include <vector>
#include <chrono>

#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::function;
using std::vector;
using std::pair;

using Eigen::VectorXd;

static vector<pair<int,double>> DEFAULT_VECTOR;

class MCMC_generator 
{
public:
	int DIM; // dimension of vector
	double alpha; // radius of multidimensional ball to jump to
	int subchain_length;
	vector<pair<int, double>> to_wrap;

	bool burnin_done;

	VectorXd current_point{ DIM };

	function<double(VectorXd)> f;

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generator{ seed }; 
	
	void set_initial_point( std::vector<double> ip );

	double nextDouble( const double& min, const double& max );
	void nextGaussianVec( VectorXd &v, VectorXd &mean );
	double wrapMax( const double& x, const double& max );

	void burnin( VectorXd initial_point, const int& burnin_length );
	VectorXd metro_step( VectorXd& x );

	VectorXd generate_point( );

	MCMC_generator( 
		function<double(VectorXd)> f, 
		const int& DIM, 
		const double& alpha, 
		const int& subchain_length, 
		vector<pair<int, double>>& to_wrap = DEFAULT_VECTOR 
			);
	~MCMC_generator();
};
