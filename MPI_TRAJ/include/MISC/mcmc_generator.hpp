#pragma once

#include <iostream>
#include <random>
#include <functional>
#include <vector>

#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::function;
using std::vector;
using std::pair;

using Eigen::VectorXf;

class MCMC_generator 
{
public:
	int DIM; // dimension of vector
	double alpha; // radius of multidimensional ball to jump to
	int n; // number of burnin steps to do before point selection 

	function<double(VectorXf)> f;
	std::mt19937 generator;

	double nextDouble( const double& min, const double& max );
	void nextGaussianVec( VectorXf &v, VectorXf &mean );
	double wrapMax( const double& x, const double& max );

	VectorXf metro_step( VectorXf& x );

	VectorXf generate_point(
		vector<pair<int, double>>& to_wrap = DEFAULT_VECTOR
						   );

	MCMC_generator( function<double(VectorXf)> f, const int& DIM, const double& alpha );
	~MCMC();
};
