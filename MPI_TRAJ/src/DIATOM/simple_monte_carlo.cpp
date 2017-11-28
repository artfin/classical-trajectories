#include <mpi.h>

#include <iostream>
#include <random>

// matrix multiplication
#include "matrix.h"

// FileReader class
#include "file.h"
// GridParameters class
#include "gridparameters.h"

// should be included BEFORE Gear header files
// due to namespace overlap
#include <vector>

// miscellaneous functions
#include "fft.h"
#include "gsl_custom.h"
#include "aux.h"

// physical constants
#include "constants.h"

// Gear header files
#include "basis.h"
#include "vmblock.h"
#include "gear.h"
#include "t_dgls.h"

// library for FFT
#include <fftw3.h>

#include <ctime>

// ############################################
const int MaxTrajectoryLength = 150000;
const double FREQ_MAX = 700.0;
// ############################################

// ############################################
// Exit tag for killing slave
const int EXIT_TAG = 42;
// ############################################

// ############################################
// Number of initial conditions per trajectory
const int ICPERTRAJ = 6; 
// R0, b, v0 -- atomic units (ALU, ALU, AVU)
// counter, R0, b, b_weight(integration), v0, v0_weight(integration)  
// ############################################

// ############################################
const double Temperature = 295; 
// ############################################

// ############################################
const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double PROTON_TO_ELECTRON_RATIO = 1836.15267389; 
const double MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 
// ############################################

// ############################################
double V0_STEP;
double B_STEP;
// ############################################

// ############################################
// sampling time in atomic time units
const double sampling_time = 100; 
// ############################################

using namespace std;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}

void show( string name, vector<double> v )
{
	cout << "#########################" << endl;
	
	for ( int i = 0; i < v.size(); i++ )
	{
		cout << name << "[" << i << "] = " << v[i] << endl;
	}
	
	cout << "#########################" << endl;
}

void show( string name1, string name2, vector<double> v1, vector<double> v2 )
{
	cout << "#########################" << endl;

	for ( int i = 0; i < v1.size(); i++ )
	{
		cout << name1 << "[" << i << "] = " << v1[i] << "; " <<
				name2 << "[" << i << "] = " << v2[i] << endl;
	}	
	
	cout << "#########################" << endl;
}

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

void save( vector<double> v1, vector<double> v2, string filename )
{
	assert( v1.size() == v2.size() && "Sizes of saved vectors should be equal" );

	ofstream file( filename );

	for ( int k = 0; k < v1.size(); k++ )
	{
		file << v1[k] << " " << v2[k] << endl;
	}
	
	file.close();	
}

void save( double m2, string filename )
{
	ofstream file( filename );

	file << "M2: " << m2 << endl;
		
	file.close();
}

vector<double> create_frequencies_vector( void )
{
	double FREQ_STEP = 1.0 / (sampling_time * constants::ATU) / constants::LIGHTSPEED_CM / MaxTrajectoryLength; // cm^-1
	//cout << "FREQ_STEP: " << FREQ_STEP << endl;
	int FREQ_SIZE = (int) FREQ_MAX / FREQ_STEP + 1;

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

	// initialzing GridParameters class
	// maybe a simple struct would do??? will see
	GridParameters grid;
	// initializing file reader class and passing it a grid class to fill in
	FileReader fileReader( "grid.in", grid ); 
}

void copy_initial_conditions( double* ics, REAL* y0, const int length )
{
	for ( int i = 0; i < length; i++ )
	{
		y0[i] = ics[i + 1];
	}
}

void slave_code( int world_rank )
{
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
	int  i;         // loop counter

	void *vmblock;  // List of dynamically allocated vectors

	// array to store initial conditions
	double *ics = new double [ICPERTRAJ];

	// auxiliary variable to store status of message
	MPI_Status status;

	// receiving integration type
  	int l;
	MPI_Recv( &l, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	//cout << "Received length of buffer: " << l << endl;

	char *buf = new char [l];
	MPI_Recv( buf, l, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
	// converting char pointer to string
	string integration_type( buf, l );
	delete [] buf;

	int integration_denumerator = 0;
	if ( integration_type == "simpson" )
	{
		//cout << "Simpson integration type. Setting denumerator to 3." << endl;
		integration_denumerator = 3;
	}
	if ( integration_type == "trapezoid" )
	{
		//cout << "Trapezoid integartion type. Setting denumerator to 2." << endl;
		integration_denumerator = 2;
	}
	if ( integration_denumerator == 0 )
	{
		cout << "UNKNOWN INTEGRATION TYPE" << endl;
		exit( 1 );
	}

	// receiving B_STEP and V0_STEP
	MPI_Recv( &B_STEP, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
	MPI_Recv( &V0_STEP, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

	//cout << "Received B_STEP: " << B_STEP << endl;
	//cout << "Received V0_STEP: " << V0_STEP << endl;
	
	// initializing special fourier class
	Fourier fourier( MaxTrajectoryLength );
	
	vector<double> freqs = create_frequencies_vector( );
	int FREQ_SIZE = freqs.size();
	double FREQ_STEP = freqs[1] - freqs[0];

	double specfunc_coeff = 1.0/(4.0*M_PI)/constants::EPSILON0 * pow(sampling_time * constants::ATU, 2)/2.0/M_PI * pow(constants::ADIPMOMU, 2);;

	double spectrum_coeff = 8.0*M_PI*M_PI*M_PI/3.0/constants::PLANCKCONST/constants::LIGHTSPEED * 1.0/4.0/M_PI/constants::EPSILON0 * pow(sampling_time * constants::ATU, 2)/2.0/M_PI * pow(constants::ADIPMOMU, 2) * pow(constants::LOSHMIDT_CONSTANT, 2);	
		
	// j -> erg; m -> cm
	double SPECFUNC_POWERS_OF_TEN = 1e19;
	// m^-1 -> cm^-1
	double SPECTRUM_POWERS_OF_TEN = 1e-2;

	while ( true )
	{
		MPI_Recv(&ics[0], ICPERTRAJ, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		//cout << ">> Process " << world_rank << " received new initial conditions." << endl; 

		if ( status.MPI_TAG == EXIT_TAG )
		{
			break;
		}

		//cout << "Received conditions: " << endl;
		//cout << "ics[0] = " << ics[0] << endl;
		//cout << "ics[1] = " << ics[1] << endl;
		//cout << "ics[2] = " << ics[2] << endl;
		//cout << "ics[3] = " << ics[3] << endl;
		//cout << "ics[4] = " << ics[4] << endl;
		//cout << "ics[5] = " << ics[5] << endl;
		//cout << "########################" << endl;

		int trajectory_number = ics[0];
		double R = ics[1];
		double b = ics[2];
		double b_weight_integration = ics[3];
		double v0 = -ics[4];
		double v0_weight_integration = ics[5];

		double pR = MU * v0;
		double pT = b * pR;
	
		// initial orientation doesn't matter 
		// it gives the same result as if theta is set to theta = 0 	
		// theta is choosed uniformly from [0, 2pi]
		double theta = nextDouble( 0.0, 2 * M_PI );	

		N = 4;
		vmblock = vminit();
		y0 = (REAL*) vmalloc(vmblock, VEKTOR, N, 0);

		y0[0] = R;
		y0[1] = pR;
		y0[2] = theta;
		y0[3] = pT;

		// out of memory?
		if ( !vmcomplete(vmblock) )
		{ 
			printf("mgear: out of memory.\n");
			return;
		}

		// according to Ivanov:
		// sampling time = 50 fs

		epsabs = 1E-13;
		epsrel = 1E-13;

		t0 = 0.0;

		h = 0.1;         		// initial step size
		xend = sampling_time;   // initial right bound of integration
		fmax = 1e8;  	 		// maximal number of calls 

		//cout << "y0[0] = " << y0[0] << endl;
		//cout << "y0[1] = " << y0[1] << endl;
		//cout << "y0[2] = " << y0[2] << endl;
		//cout << "y0[3] = " << y0[3] << endl;
	
		double energy = pow(y0[1], 2) / 2 / MU + pow(y0[3], 2) / ( 2 * MU * pow(y0[0], 2));
		//cout << "energy: " << energy << endl;
		
		cout << "b (a.u): " << b << "; b (SI): " << b * constants::ALU << endl;
		cout << "v0 (a.u.): " << v0 << "; v0 (SI): " << v0 * constants::AVU << endl;

		double b_weight = fabs( 2 * M_PI * v0 * b ) * constants::AVU * constants::ALU * b_weight_integration;
		cout << "b_weight (SI): " << b_weight << endl;
		// integration-denumerator = 2.0 -- trapezoid rule
		// integration-denumerator = 3.0 -- simpson rule

		double v0_weight_constant = pow((MU * constants::AMU / (2 * M_PI * constants::BOLTZCONST * Temperature)), 1.5);
		double v0_weight = v0_weight_constant * exp(- MU * constants::AMU * pow(v0 * constants::AVU, 2) / ( 2 * constants::BOLTZCONST * Temperature )) * 4 * M_PI * pow(v0, 2) * V0_STEP * pow(constants::AVU, 3) * v0_weight_integration;
		cout << "v0_weight: " << v0_weight << endl;

		//double fourier_weight = pow(sampling_time * constants::ATU, 2);
		double stat_weight = b_weight * v0_weight * B_STEP * constants::ALU;
		cout << "stat_weight: " << stat_weight << endl;

		double stat_weight_si = stat_weight * constants::ALU * constants::ALU * constants::AVU;
		cout << "stat_weight_si: " << stat_weight_si << endl;

		//cout << "B_STEP: " << B_STEP << endl;
		//cout << "V0_STEP: " << V0_STEP << endl;
		//cout << "b_weight: " << b_weight << endl;
		//cout << "v0_weight: " << v0_weight << endl;
		//cout << "weight: " << weight << endl;

		int counter = 0;
		double end_value = y0[0] + 0.001; 

		// dipole moment in laboratory frame
		vector<double> temp( 3 );

		vector<double> dipx;
		vector<double> dipy;
		vector<double> dipz;

		while ( y0[0] < end_value ) 
		{
			fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);

			if ( fehler != 0 ) 
			{
				printf(" Gear4: error n° %d\n", 10 + fehler);
				break;
			}

			//cout << "t0: " << t0 << "; r: " << y0[0] << endl;
			
			// y0[0] -- R
			// y0[1] -- PR
			// y0[2] -- \theta
			// y0[3] -- p_\theta
			transform_dipole( temp, y0[0], y0[2] );

			dipx.push_back( temp[0] );
			dipy.push_back( temp[1] );
			dipz.push_back( temp[2] );

			xend = sampling_time * (counter + 2);
			aufrufe = 0;  // actual number of calls

			counter++;
		}
		
		// length of dipole vector = number of samples
		int npoints = dipz.size();

		vector<double> specfunc;
		vector<double> spectrum;

		double specfunc_value;
		double spectrum_value;

		double uniform_integrated_spectrum = 0;
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

			specfunc_value = SPECFUNC_POWERS_OF_TEN * specfunc_coeff * stat_weight * dipfft;
			specfunc.push_back( specfunc_value );	
				
			spectrum_value = SPECTRUM_POWERS_OF_TEN * spectrum_coeff * stat_weight * omega * ( 1.0 - exp( - constants::PLANCKCONST_REDUCED * omega / ( constants::BOLTZCONST * Temperature )) ) * dipfft;
				
			uniform_integrated_spectrum += spectrum_value * FREQ_STEP; 
			// PLANCKCONST/2/PI = PLANCKCONST_REDUCED

			spectrum.push_back( spectrum_value );
		}

		cout << ">> Processing " << trajectory_number << " trajectory. npoints = " << npoints << endl;
			
		//// #################################################
		// Sending data
		MPI_Send( &specfunc[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		MPI_Send( &spectrum[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		MPI_Send( &uniform_integrated_spectrum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );	
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
