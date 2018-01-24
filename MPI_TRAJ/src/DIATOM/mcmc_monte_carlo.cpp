#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include <functional>

#include <Eigen/Dense>

// matrix multiplication
#include "matrix.h"

// FileReader class
#include "file.h"
// Parameters class
#include "parameters.h"
// SpectrumInfo class
#include "spectrum_info.hpp"
// MCMC_generator class
#include "mcmc_generator.hpp"
// Integrand class
#include "integrand.hpp"
// Integrator class
#include "integrator.hpp"

// should be included BEFORE Gear header files
// due to namespace overlap
#include <vector>

// miscellaneous functions
#include "fft.h"

// physical constants
#include "constants.h"

// library for FFT
#include <fftw3.h>

// Gear header files
#include "basis.h"
#include "vmblock.h"
#include "gear.h"
#include "t_dgls.h"

using namespace std;
using namespace std::placeholders;
using Eigen::VectorXd;

// ############################################
// Exit tag for killing slave
const int EXIT_TAG = 42;
// ############################################

// ############################################
const double DALTON_UNIT = 1.660539040 * 1e-27;

const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double MU_SI = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * DALTON_UNIT; 
const double MU = MU_SI / constants::AMU;

static const double MYPI = 3.141592653589793; 
const double TWO_PI = 2 * MYPI; 
// ############################################

void syst (REAL t, REAL *y, REAL *f)
{
  	(void)(t); // avoid unused parameter warning 

	double *out = new double[4];

	rhs( out, y[0], y[1], y[2], y[3] );

	//cout << "out[0]: " << out[0] << endl;
	//cout << "out[1]: " << out[1] << endl;
	//cout << "out[2]: " << out[2] << endl;
	//cout << "out[3]: " << out[3] << endl;

	f[0] = out[0]; // \dot{R} 
	f[1] = out[1]; // \dot{p_R}
	f[2] = out[2]; // \dot{\theta}
	f[3] = out[3]; // \dot{p_\theta}

	delete [] out;
}

// x = [ pR, theta, pT ]
double target_distribution( VectorXd x, const double& R, const double& temperature )
{
	double pR = x( 0 );
	double theta = x( 1 );
	double pT = x( 2 );

	double h = pow(pR, 2) / (2 * MU) + pow(pT, 2) / (2 * MU * R * R);
	return exp( -h * constants::HTOJ / (constants::BOLTZCONST * temperature )); 
}

vector<double> create_frequencies_vector( Parameters& parameters )
{
	double FREQ_STEP = 1.0 / (parameters.sampling_time * constants::ATU) / constants::LIGHTSPEED_CM / parameters.MaxTrajectoryLength; // cm^-1
	//cout << "FREQ_STEP: " << FREQ_STEP << endl;
	int FREQ_SIZE = (int) parameters.FREQ_MAX / FREQ_STEP + 1;

	vector<double> freqs( FREQ_SIZE );
	for(int k = 0; k <  FREQ_SIZE; k++) 
	{
		freqs[k] = k * FREQ_STEP;
	}

	return freqs;
}

void master_code( int world_size )
{
	MPI_Status status;
	int source;

	Parameters parameters;
	FileReader fileReader( "parameters.in", &parameters ); 
	//parameters.show_parameters();

	function<double(VectorXd)> target = bind( target_distribution, _1, parameters.RDIST, parameters.Temperature );
	
	Integrand integrand( target, 3 );
	integrand.set_limits()->add_limit( 0, "-inf", "+inf" )
						  ->add_limit( 1, 0.0, 2.0 * M_PI )
						  ->add_limit( 2, "-inf", "+inf" );
	
	Integrator integrator( integrand, 0 );
	integrator.set_callback();
	
	int niter = 10, ndots = 50000;	
	double ham_integral = integrator.run_integration( niter, ndots );
	ham_integral = ham_integral * constants::HTOJ * constants::HTOJ * constants::ATU / constants::ALU;	

	int sent = 0;	
	int received = 0;

	// status of calculation
	bool is_finished = false;

	vector<double> freqs = create_frequencies_vector( parameters );
	int FREQ_SIZE = freqs.size();

	// creating objects to hold spectrum info
	SpectrumInfo classical( FREQ_SIZE, "classical" );
	SpectrumInfo d1( FREQ_SIZE, "d1" );
	SpectrumInfo d2( FREQ_SIZE, "d2" );
	SpectrumInfo d3( FREQ_SIZE, "d3" );
	
	// wrapping second argument (argument 1):
	pair<int, double> p1(1, 2*M_PI); 
	vector<pair<int, double>> to_wrap{ p1 };

	MCMC_generator generator( target, parameters.DIM, parameters.alpha, parameters.subchain_length, to_wrap );
	
	VectorXd initial_point = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( parameters.initial_point.data(), parameters.initial_point.size());	
	generator.burnin( initial_point, 100000 );

	VectorXd p;
	// sending first trajectory	
	for ( int i = 1; i < world_size; i++ )
	{
		p = generator.generate_point();
		MPI_Send( p.data(), parameters.DIM, MPI_DOUBLE, i, 0, MPI_COMM_WORLD );
		MPI_Send( &sent, 1, MPI_INT, i, 0, MPI_COMM_WORLD );

		sent++;
	}

	clock_t start = clock();

	while( true )
	{
		if ( is_finished )	
		{
			for ( int i = 1; i < world_size; i++ )
			{	
				MPI_Send( &is_finished, 1, MPI_INT, i, EXIT_TAG, MPI_COMM_WORLD );
			}

			break;
		}

		// ############################################################
		// Receiving data
		classical.receive( source, true );
		received++;
		// ############################################################

		classical.add_package_to_total();

		string name = "temp";
		stringstream ss;
		if ( received % 10 == 0 )
		{
			ss << received;
			classical.saving_procedure( parameters, freqs, name + ss.str() + ".txt", "total" );
			classical.zero_out_total();
		}

		if ( received == parameters.NPOINTS )
		{
			double multiplier = 1.0 / parameters.NPOINTS; 
			//cout << "multiplier: " << multiplier << endl;

			classical.multiply_total( multiplier / ham_integral );

			cout << ">>Saving spectrum" << endl << endl;
			classical.saving_procedure( parameters, freqs ); 

			is_finished = true;
		}		

		if ( sent < parameters.NPOINTS )
		{
			p = generator.generate_point();
			MPI_Send( p.data(), parameters.DIM, MPI_DOUBLE, source, 0, MPI_COMM_WORLD );
			MPI_Send( &sent, 1, MPI_INT, source, 0, MPI_COMM_WORLD );

			sent++;
		}
	}
}

//double d1_corrector( double omega, double kT )
//{
	//return 2.0 / (1.0 + exp(-constants::PLANCKCONST_REDUCED * omega / kT));
//}	

//double d2_corrector( double omega, double kT )
//{
	//return constants::PLANCKCONST_REDUCED * omega / kT / (1.0 - exp(-constants::PLANCKCONST_REDUCED * omega / kT));
//}

//double d3_corrector( double omega, double kT )
//{
	//return exp( constants::PLANCKCONST_REDUCED * omega / kT / 2.0 );
//}

void slave_code( int world_rank )
{
	// it would be easier for slave to read parameters file by himself, rather than sending him parameters object (or just several parameters)
	Parameters parameters;
	FileReader fileReader( "parameters.in", &parameters ); 

	MPI_Status status;
	
	// #####################################################
	// Integrating hamiltonian over the phase space
	function<double(VectorXd)> target = bind( target_distribution, _1, parameters.RDIST, parameters.Temperature );
	
	Integrand integrand( target, 3 );
	integrand.set_limits()->add_limit( 0, "-inf", "+inf" )
						  ->add_limit( 1, 0.0, 2.0 * M_PI )
						  ->add_limit( 2, "-inf", "+inf" );
	
	Integrator integrator( integrand, world_rank );
	integrator.set_callback();
	
	int niter = 10, ndots = 50000;	
	double ham_integral = integrator.run_integration( niter, ndots );
	// #####################################################

	// #####################################################
	// initializing special fourier class
	Fourier fourier( parameters.MaxTrajectoryLength );

	vector<double> freqs = create_frequencies_vector( parameters );
	int FREQ_SIZE = freqs.size();
	double FREQ_STEP = freqs[1] - freqs[0];

	double specfunc_coeff = 1.0/(4.0*M_PI)/constants::EPSILON0 * pow(parameters.sampling_time * constants::ATU, 2)/2.0/M_PI * pow(constants::ADIPMOMU, 2);

	double spectrum_coeff = 8.0*M_PI*M_PI*M_PI/3.0/constants::PLANCKCONST/constants::LIGHTSPEED * 1.0/4.0/M_PI/constants::EPSILON0 * pow(parameters.sampling_time * constants::ATU, 2)/2.0/M_PI * pow(constants::ADIPMOMU, 2) * pow(constants::LOSHMIDT_CONSTANT, 2);	

	// j -> erg; m -> cm
	double SPECFUNC_POWERS_OF_TEN = 1e19;
	// m^-1 -> cm^-1
	double SPECTRUM_POWERS_OF_TEN = 1e-2;

	double kT = constants::BOLTZCONST * parameters.Temperature;
	// #####################################################
	
	//// #####################################################
	REAL epsabs;    //  absolute error bound
	REAL epsrel;    //  relative error bound    
	REAL t0;        // left edge of integration interval
	REAL *y0;       // [0..n-1]-vector: initial value, approxim. 
	REAL h;         // initial, final step size
	REAL xend;      // right edge of integration interval 

	long fmax;      // maximal number of calls of right side in gear4()
	long aufrufe;   // actual number of function calls
	int  N;         // number of DEs in system
	int  fehler;    // error code from umleiten(), gear4()

	void *vmblock;  // List of dynamically allocated vectors


	N = 4;
	vmblock = vminit();
	y0 = (REAL*) vmalloc(vmblock, VEKTOR, N, 0);

	// accuracy of trajectory
	epsabs = 1E-13;
	epsrel = 1E-13;

	fmax = 1e8;  		  // maximal number of calls 
	// #####################################################

	// creating objects to hold spectal info
	SpectrumInfo classical;

	vector<double> p( parameters.DIM );
	int traj_counter = 0;

	while ( true )
	{
		MPI_Recv( &p[0], parameters.DIM, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
		
		if ( status.MPI_TAG == EXIT_TAG )
		{
			cout << "Received exit tag." << endl;
			break;
		}
		MPI_Recv( &traj_counter, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

		y0[0] = parameters.RDIST; 
		y0[1] = p[0];
		y0[2] = p[1]; 
		y0[3] = p[2]; 

		//cout << "#####" << endl;
		//cout << "ics.R = " << y0[0] << endl;
		//cout << "ics.pR = " << y0[1] << endl;
		//cout << "ics.theta = " << y0[2] << endl;
		//cout << "ics.pT = " << y0[3] << endl;
		//cout << "#####" << endl;

		// out of memory?
		//if ( !vmcomplete(vmblock) )
		//{ 
		//cout << "mgear: out of memory" << endl;
		//return;
		//}

		int counter = 0;
		double R_end_value = y0[0] + 0.001;

		// dipole moment in laboratory frame
		vector<double> temp( 3 );
		vector<double> dipx;
		vector<double> dipy;
		vector<double> dipz;

		// #####################################################
		t0 = 0.0;

		h = 0.1;    		  // initial step size
		xend = parameters.sampling_time; // initial right bound of integration
		// #####################################################

		clock_t start = clock();

		// #####################################################
		while ( y0[0] < R_end_value ) 
		{
			if ( counter == parameters.MaxTrajectoryLength )
			{
				//cout << "Trajectory cut!" << endl;
				break;
			}

			fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);

			//cout << "%%%" << endl;
			//cout << "t0: " << t0 << endl;
			//cout << "xend: " << xend << endl;
			//cout << "%%%" << endl;

			if ( fehler != 0 ) 
			{
				cout << "Gear4: error n = " << 10 + fehler << endl;
				break;
			}

			//cout << "t0: " << t0 << "; r: " << y0[0] << endl;

			transform_dipole( temp, y0[0], y0[2] );

			//cout << "t: " << t0 * constants::ATU <<
				//"; R (alu): " << y0[0] << 	
				//"; R: " << y0[0] * constants::ALU << 
				//"; dipole z: " << temp[2] * constants::ADIPMOMU << endl;

			dipx.push_back( temp[0] );
			dipy.push_back( temp[1] );
			dipz.push_back( temp[2] );

			xend = parameters.sampling_time * (counter + 2);

			aufrufe = 0;  // actual number of calls

			counter++;
		}
		// #####################################################

		// #####################################################
		// length of dipole vector = number of samples
		int npoints = dipz.size();

		//cout << "FREQ_SIZE: " << FREQ_SIZE << endl;

		// zeroing input arrays
		fourier.zero_out_input( );

		// copying data to them
		copy_to( dipx, fourier.inx );
		copy_to( dipy, fourier.iny );
		copy_to( dipz, fourier.inz );

		// executing fourier transform
		fourier.do_fourier( );

		double omega, dipfft;

		double ReFx, ReFy, ReFz;
		double ImFx, ImFy, ImFz;

		double specfunc_value_classical;
		double spectrum_value_classical;

		for ( int k = 0; k < FREQ_SIZE; k++ )
		{
			omega = 2.0 * M_PI * constants::LIGHTSPEED_CM * freqs[k];

			ReFx = fourier.outx[k][0];
			ReFy = fourier.outy[k][0];
			ReFz = fourier.outz[k][0];
			ImFx = fourier.outx[k][1];
			ImFy = fourier.outy[k][1];
			ImFz = fourier.outz[k][1];

			dipfft = ReFx * ReFx + ReFy * ReFy + ReFz * ReFz +
 					 ImFx * ImFx + ImFy * ImFy + ImFz * ImFz;
			//cout << "dipfft[" << k << "] = " << dipfft * constants::ADIPMOMU * constants::ADIPMOMU << endl; 

			specfunc_value_classical = SPECFUNC_POWERS_OF_TEN * specfunc_coeff * dipfft;
			classical.specfunc_package.push_back( specfunc_value_classical );

			spectrum_value_classical = SPECTRUM_POWERS_OF_TEN * spectrum_coeff * omega *  ( 1.0 - exp( - constants::PLANCKCONST_REDUCED * omega / kT ) ) * dipfft;
			classical.spectrum_package.push_back( spectrum_value_classical );

			classical.m2_package += spectrum_value_classical * FREQ_STEP; 
		}

		cout << "(" << world_rank << ") Processing " << traj_counter << " trajectory. npoints = " << npoints << "; time = " << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;

		//// #################################################
		// Sending data
		classical.send();
		classical.clear_package();
		// #################################################
	}
}

int main( int argc, char* argv[] )
{
	//Initialize the MPI environment
	MPI_Init( &argc, &argv );

	//getting id of the current process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	//getting number of running processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); 

	if ( world_rank == 0 ) 
	{
		clock_t start = clock();

		master_code( world_size );

		cout << "Time elapsed: " << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
	}	
	else
	{
		slave_code( world_rank );
	}

	MPI_Finalize();

	return 0;
}
