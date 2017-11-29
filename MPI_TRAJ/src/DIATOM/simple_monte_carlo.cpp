#include <mpi.h>

#include <iostream>
#include <random>

// matrix multiplication
#include "matrix.h"

// FileReader class
#include "file.h"
// Parameters class
#include "parameters.h"

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
// initial R (A.U.)
const double RDIST = 40.0;
// ############################################

// ############################################
const double Temperature = 295; 
// ############################################

// ############################################
const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double PROTON_TO_ELECTRON_RATIO = 1836.15267389; 
const double MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 
const double MU_SI = MU * constants::AMU;

const double TWO_PI = 2 * M_PI;
// ############################################

// ############################################
double V0_STEP;
double B_STEP;
// ############################################

// ############################################
// sampling time in atomic time units
const double sampling_time = 100; 
// ############################################

struct ICHamPoint
{
	double R;
	double pR;
	double theta;
	double pT;
};

using namespace std;

static mt19937 uniform_generator;

static double UniformDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( uniform_generator );
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

	Parameters parameters;

	FileReader fileReader( "parameters.in", &parameters ); 
	parameters.show_parameters();

	// #####################################################
	// creating custom data type for MPI
	int blocksCount = 3; // number of entities inside struct
	int blocksLength[3] = {1, 1, 1}; // lengths of entities inside struct

	// types of variables in struct
	MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};

	MPI_Aint offsets[3]; // mpi_aint holds address (pointer)
	MPI_Datatype MPI_ICPoint;
	offsets[0] = offsetof( ICPoint, counter ); // sizes of variables?
	offsets[1] = offsetof( ICPoint, b );
	offsets[2] = offsetof( ICPoint, v0 );

	MPI_Type_create_struct( blocksCount, blocksLength, offsets, types, &MPI_ICPoint );
	MPI_Type_commit( &MPI_ICPoint );
	// #####################################################

	ICPoint p;
	int point_counter = 0;
	int cycle_counter = 0;

	// sending first point to slaves
	for ( int i = 1; i < world_size; i++ )
	{
		p = parameters.generate_uniform_point( );
		p.counter = point_counter;
			
		MPI_Send( &p, 1, MPI_ICPoint, i, 0, MPI_COMM_WORLD );
		//cout << "Master sends first message" << endl;

		point_counter++;
	}

	// status of calculation: set to true if all cycles are done
	bool finished = false;

	// number of alive processes
	int alive = world_size - 1;
	
	vector<double> freqs = create_frequencies_vector( );
	int FREQ_SIZE = freqs.size();
	
	vector<double> specfunc_package( FREQ_SIZE );
	vector<double> total_specfunc( FREQ_SIZE );
	vector<double> spectrum_package( FREQ_SIZE );
	vector<double> total_spectrum( FREQ_SIZE );
	
	double uniform_integrated_spectrum_package;
	double uniform_integrated_spectrum_total = 0.0;

	while( true )
	{
		if ( alive == 0 )
		{
			//cout << "All slaves are dead" << endl;
			break;
		}

		if ( point_counter == parameters.CYCLE_POINTS )
		{
			cout << "Cycle " << cycle_counter << " is finished." << endl;
	
			// multiplyting spectrum and spectral function by AREA and dividing by NPOINTS
			double b_length = parameters.B_MAX - parameters.B_MIN;
			double v0_length = parameters.V0_MAX - parameters.V0_MIN;
			double multiplier = b_length * v0_length / point_counter;
			for ( int i = 0; i < total_specfunc.size(); i++ )
			{
				total_specfunc[i] *= multiplier; 
				total_spectrum[i] *= multiplier; 
			}

			uniform_integrated_spectrum_total *= multiplier; 

			cout << "Saving data to files..." << endl;
					
			save( freqs, total_specfunc, parameters.specfunc_filename );
			save( freqs, total_spectrum, parameters.spectrum_filename );
			save( uniform_integrated_spectrum_total, parameters.m2_filename );

			point_counter = 0;
			cycle_counter++;
		}		
		
		if ( cycle_counter == parameters.CYCLES )
		{
			finished = true;
		}
	
		// ############################################################
		// Receiving data
		MPI_Recv( &specfunc_package[0], FREQ_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		source = status.MPI_SOURCE;

		MPI_Recv( &spectrum_package[0], FREQ_SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		MPI_Recv (&uniform_integrated_spectrum_package, 1, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		// ############################################################
		
		// ############################################################
		for ( int i = 0; i < FREQ_SIZE; i++ )
		{
			total_specfunc[i] += specfunc_package[i];
			total_spectrum[i] += spectrum_package[i];
		}
		uniform_integrated_spectrum_total += uniform_integrated_spectrum_package;
		// ############################################################
		
		if ( !finished )
		{
			p = parameters.generate_uniform_point();
			p.counter = point_counter;	
			point_counter++;
			
			MPI_Send( &p, 1, MPI_ICPoint, source, 0, MPI_COMM_WORLD );
		}
		else
		{
			MPI_Send( &finished, 1, MPI_INT, source, EXIT_TAG, MPI_COMM_WORLD );
			alive--;
		}
	}
}

void copy_initial_conditions( double* ics, REAL* y0, const int length )
{
	for ( int i = 0; i < length; i++ )
	{
		y0[i] = ics[i + 1];
	}
}

void show_point( ICPoint p )
{
	cout << "%%%%%%%%%%" << endl;
	cout << "p.counter: " << p.counter << endl;
	cout << "p.b: " << p.b << endl;
	cout << "p.v0: " << p.v0 << endl;
	cout << "%%%%%%%%%%" << endl;
}

void transform_ICPoint_to_ICHamPoint( ICPoint &p, ICHamPoint &ics )
{
	// input: p.b, p.vo -- in SI 
	// output: ics.R, pR, theta, pT -- in A.U.
	
	ics.R = RDIST;
   	ics.pR = MU * p.v0 / constants::AVU;
	ics.theta = UniformDouble( 0.0, 2 * M_PI );
	ics.pT = p.b / constants::ALU * ics.pR;	
	
	//cout << "p.b: " << p.b << endl;
	//cout << "ics.pR: " << ics.pR << endl;
	//cout << "ics.pT: " << ics.pT << endl;
}

void slave_code( int world_rank )
{
	MPI_Status status;
	
	// #####################################################
	// creating custom data type for MPI
	int blocksCount = 3; // number of entities inside struct
	int blocksLength[3] = {1, 1, 1}; // lengths of entities inside struct

	// types of variables in struct
	MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};

	MPI_Aint offsets[3]; // mpi_aint holds address (pointer)
	MPI_Datatype MPI_ICPoint;
	offsets[0] = offsetof( ICPoint, counter ); // sizes of variables?
	offsets[1] = offsetof( ICPoint, b );
	offsets[2] = offsetof( ICPoint, v0 );

	MPI_Type_create_struct( blocksCount, blocksLength, offsets, types, &MPI_ICPoint );
	MPI_Type_commit( &MPI_ICPoint );
	// #####################################################
	
	// #####################################################
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
	// #####################################################
	
	// #####################################################
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

	ICPoint p;
	ICHamPoint ics;

	while ( true )
	{
		MPI_Recv( &p, 1, MPI_ICPoint, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
		//cout << "Received ICPoint: " << p.counter << endl;
		//show_point( p );
		transform_ICPoint_to_ICHamPoint( p, ics );

		if ( status.MPI_TAG == EXIT_TAG )
		{
			cout << "Received exit tag." << endl;
			break;
		}

		y0[0] = ics.R;
		y0[1] = - ics.pR;
		y0[2] = ics.theta;
		y0[3] = ics.pT;
		//cout << "y0[0] = " << y0[0] << endl;
		//cout << "y0[1] = " << y0[1] << endl;
		//cout << "y0[2] = " << y0[2] << endl;
		//cout << "y0[3] = " << y0[3] << endl;

		// out of memory?
		if ( !vmcomplete(vmblock) )
		{ 
			cout << "mgear: out of memory" << endl;
			return;
		}

		// #####################################################
		// p.v0, p.b -- in SI
		double b_integrand = 2 * M_PI * p.v0 * p.b;  
		double v0_integrand = pow( MU_SI / (TWO_PI * constants::BOLTZCONST * Temperature), 1.5 ) * exp( - MU_SI * p.v0 * p.v0 / (2 * constants::BOLTZCONST * Temperature)) * 4 * M_PI * p.v0 * p.v0; 
		
		double stat_weight = b_integrand * v0_integrand;
		// #####################################################
		
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
		xend = sampling_time; // initial right bound of integration
		// #####################################################

		// #####################################################
		while ( y0[0] < R_end_value ) 
		{
			fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);
			//cout << "%%%" << endl;
			//cout << "xend: " << xend << endl;
			//cout << "%%%" << endl;

			if ( fehler != 0 ) 
			{
				cout << "Gear4: error n = " << 10 + fehler << endl;
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
		// #####################################################
		
		// #####################################################
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

		cout << ">> Processing " << p.counter << " trajectory. npoints = " << npoints << endl;
		
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
