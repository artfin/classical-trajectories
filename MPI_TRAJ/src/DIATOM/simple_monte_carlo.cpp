#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>

// FileReader class
#include "file.h"
// Parameters class
#include "parameters.h"
// SpectrumInfo class
#include "spectrum_info.hpp"

// should be included BEFORE Gear header files
// due to namespace overlap
#include <vector>

// miscellaneous functions
#include "fft.h"

// physical constants
#include "constants.h"

// library for FFT
#include <fftw3.h>

// Trajectory class (INCLUDED LAST!)
#include "trajectory.hpp"

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

struct ICHamPoint
{
	double R;
	double pR;
	double theta;
	double pT;
};

using namespace std;

static mt19937 uniform_generator( 27717 );

static double UniformDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( uniform_generator );
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

void create_chunks( Parameters& parameters, vector<double>& B_CHUNKS_VECTOR, vector<double>& V0_CHUNKS_VECTOR )
{
	double B_CHUNK, V0_CHUNK;

	if ( parameters.B_PARTS != 1 )
	{
		B_CHUNK = (parameters.B_MAX - parameters.B_MIN) / (parameters.B_PARTS);
	}
	else 
	{
		B_CHUNK = parameters.B_MAX - parameters.B_MIN;
	}

	if ( parameters.V0_PARTS != 1 )
	{
		V0_CHUNK = (parameters.V0_MAX - parameters.V0_MIN) / (parameters.V0_PARTS );
	}
	else
	{
		V0_CHUNK = parameters.V0_MAX - parameters.V0_MIN;
	}

	double temp;
	for ( int i = 0; i < parameters.B_PARTS + 1; i++ )
	{
		temp = parameters.B_MIN + i * B_CHUNK; 
		B_CHUNKS_VECTOR.push_back( temp );
	}	
	for ( int i = 0; i < parameters.V0_PARTS + 1; i++ )
	{
		temp = parameters.V0_MIN + i * V0_CHUNK;
		V0_CHUNKS_VECTOR.push_back( temp );
	}
}

void master_code( int world_size )
{
	MPI_Status status;
	int source;

	Parameters parameters;
	FileReader fileReader( "parameters.in", &parameters ); 
	//parameters.show_parameters();
	
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

	// ##########################################################
	int b_chunk_counter = 0;
	int v0_chunk_counter = 0;
	int b_chunk_max = parameters.B_PARTS;
	int v0_chunk_max = parameters.V0_PARTS;

	vector<double> B_CHUNKS_VECTOR;
	vector<double> V0_CHUNKS_VECTOR;
	create_chunks( parameters, B_CHUNKS_VECTOR, V0_CHUNKS_VECTOR );

	for ( int i = 0; i < B_CHUNKS_VECTOR.size(); i++ )
	{
		cout << "B_CHUNKS_VECTOR[" << i << "] = " << B_CHUNKS_VECTOR[i] << endl;
	}

	double B_CHUNK;
	double V0_CHUNK;
	// ##########################################################


	// sending first trajectory	
	for ( int i = 1; i < world_size; i++ )
	{
		p = parameters.generate_uniform_point( B_CHUNKS_VECTOR[ b_chunk_counter ],
												   B_CHUNKS_VECTOR[ b_chunk_counter + 1 ],
												   V0_CHUNKS_VECTOR[ v0_chunk_counter ],
												   V0_CHUNKS_VECTOR[ v0_chunk_counter + 1 ]);
		p.counter = sent; 
			
		//cout << "###" << endl;
		//cout << "generated p.b: " << p.b << endl;
		//cout << "generated p.v0: " << p.v0 << endl;
		//cout << "p.counter: " << p.counter << endl;
		//cout << "###" << endl;

		MPI_Send( &p, 1, MPI_ICPoint, i, 0, MPI_COMM_WORLD );
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
		d1.receive( source, false );
		d2.receive( source, false );
		d3.receive( source, false );
		received++;

		// ############################################################

		//cout << "sent: " << sent << "; received: " << received << "; NPOINTS: " << parameters.NPOINTS << endl;	
		if ( received == parameters.NPOINTS )
		{
			// multiplyting spectrum and spectral function by AREA and dividing by NPOINTS
			B_CHUNK = B_CHUNKS_VECTOR[b_chunk_counter + 1] -
					  B_CHUNKS_VECTOR[b_chunk_counter];
			V0_CHUNK = V0_CHUNKS_VECTOR[v0_chunk_counter + 1] -
					   V0_CHUNKS_VECTOR[v0_chunk_counter];
			
			double multiplier = B_CHUNK * V0_CHUNK / parameters.NPOINTS; 
			//cout << "multiplier: " << multiplier << endl;

			// ##################################################################################

			classical.multiply_chunk( multiplier );
			d1.multiply_chunk( multiplier );
			d2.multiply_chunk( multiplier );
			d3.multiply_chunk( multiplier );

			// ##################################################################################

			cout << "#####################################" << endl;
			cout << "Current chunk: (b) " << b_chunk_counter << " (v0) " << v0_chunk_counter << endl;
			cout << "B: " << B_CHUNKS_VECTOR[b_chunk_counter] << " -- " << B_CHUNKS_VECTOR[b_chunk_counter + 1] << endl;
			cout << "V0: " << V0_CHUNKS_VECTOR[v0_chunk_counter] << " -- " << V0_CHUNKS_VECTOR[v0_chunk_counter + 1] << endl;
			cout << "M2 (class): " << classical.m2_chunk << endl;
			cout << "M2 (d1): " << d1.m2_chunk << endl;
			cout << "M2 (d2): " << d2.m2_chunk << endl;
			cout << "M2 (d3): " << d3.m2_chunk << endl;
			cout << "Time for chunk: " << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
			start = clock();
			cout << "#####################################" << endl;
			
			// adding chunk spectrum to total; zeroing out chunk data
			//cout << ">> Added chunk spectrum to total" << endl << endl;	
			classical.add_chunk_to_total();
			d1.add_chunk_to_total();
			d2.add_chunk_to_total();
			d3.add_chunk_to_total();

			//cout << "point counter is set to 0" << endl << endl;
			sent = 0;
			received = 0;
			// ##################################################

			cout << ">>Saving spectrum" << endl << endl;
			classical.saving_procedure( parameters, freqs ); 
			d1.saving_procedure( parameters, freqs );
			d2.saving_procedure( parameters, freqs );
			d3.saving_procedure( parameters, freqs );

			// ##################################################
			if ( b_chunk_counter < b_chunk_max - 1 )
			{
				//cout << ">> b_chunk_counter is not maximum. incrementing..." << endl;
				b_chunk_counter++;
			}	
			else
			{
				if ( v0_chunk_counter < v0_chunk_max - 1 )
				{
					//cout << ">> v0_chunk_counter is not maximum. Incrementing..." << endl;
					//cout << ">> setting b_chunk_counter to 0" << endl;
					v0_chunk_counter++;
					b_chunk_counter = 0;
				}
				else
				{
					cout << ">> all chunks are integrated!" << endl;
					// setting finished flag to true
					is_finished = true;
					// starting new iteration which checks the finished flag 
					continue;
				}
			}
			// ##################################################

			// ##################################################
			for ( int i = 1; i < world_size; i++ )
			{
				p = parameters.generate_uniform_point( B_CHUNKS_VECTOR[ b_chunk_counter ],
												   	   B_CHUNKS_VECTOR[ b_chunk_counter + 1 ],
												   	   V0_CHUNKS_VECTOR[ v0_chunk_counter ],
													   V0_CHUNKS_VECTOR[ v0_chunk_counter + 1 ]);
				p.counter = sent; 	
			
				//cout << "###" << endl;
				//cout << "generated p.b: " << p.b << endl;
				//cout << "generated p.v0: " << p.v0 << endl;
				//cout << "p.counter: " << p.counter << endl;
				//cout << "###" << endl;

				MPI_Send( &p, 1, MPI_ICPoint, i, 0, MPI_COMM_WORLD );
				//cout << "point_counter: " << point_counter << endl;
				sent++;
			}
			// ##################################################
		}		
		
		// ############################################################
		classical.add_package_to_chunk();

		d1.add_package_to_chunk();	
		d2.add_package_to_chunk();	
		d3.add_package_to_chunk();

		//cout << "Added package spectrum to chunk." << endl;
		// ############################################################
		
		
		if ( sent < parameters.NPOINTS )
		{
			p = parameters.generate_uniform_point( B_CHUNKS_VECTOR[ b_chunk_counter ],
												   B_CHUNKS_VECTOR[ b_chunk_counter + 1 ],
												   V0_CHUNKS_VECTOR[ v0_chunk_counter ],
												   V0_CHUNKS_VECTOR[ v0_chunk_counter + 1 ]);
			p.counter = sent; 
			
			//cout << "###" << endl;
			//cout << "generated p.b: " << p.b << endl;
			//cout << "generated p.v0: " << p.v0 << endl;
			//cout << "p.counter: " << p.counter << endl;
			//cout << "###" << endl;
			
			MPI_Send( &p, 1, MPI_ICPoint, source, 0, MPI_COMM_WORLD );
			sent++;
		}
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

void transform_ICPoint_to_ICHamPoint( Parameters& parameters, ICPoint& p, ICHamPoint& ics )
{
	// input: p.b, p.vo -- in SI 
	// output: ics.R, pR, theta, pT -- in A.U.
	
	ics.R = parameters.RDIST;
   	ics.pR = MU * p.v0 / constants::AVU;
	ics.theta = UniformDouble( 0.0, 2 * M_PI );
	ics.pT = p.b / constants::ALU * ics.pR;	
	
	//cout << "p.b: " << p.b << endl;
	//cout << "ics.pR: " << ics.pR << endl;
	//cout << "ics.pT: " << ics.pT << endl;
}

double d1_corrector( double omega, double kT )
{
	return 2.0 / (1.0 + exp(-constants::PLANCKCONST_REDUCED * omega / kT));
}	

double d2_corrector( double omega, double kT )
{
	return constants::PLANCKCONST_REDUCED * omega / kT / (1.0 - exp(-constants::PLANCKCONST_REDUCED * omega / kT));
}

double d3_corrector( double omega, double kT )
{
	return exp( constants::PLANCKCONST_REDUCED * omega / kT / 2.0 );
}

void slave_code( int world_rank )
{
	// it would be easier for slave to read parameters file by himself, rather than sending him parameters object (or just several parameters)
	Parameters parameters;
	FileReader fileReader( "parameters.in", &parameters ); 

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

	ICPoint p;
	ICHamPoint ics;

	Trejectory trajectory( 4 ); // 4 equations

	// creating objects to hold spectal info
	SpectrumInfo classical;
	SpectrumInfo d1( d1_corrector );
	SpectrumInfo d2( d2_corrector );
	SpectrumInfo d3( d3_corrector );

	while ( true )
	{
		MPI_Recv( &p, 1, MPI_ICPoint, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
		//cout << "(" << world_rank << ") Received ICPoint: " << p.counter << endl;
		//show_point( p );
		transform_ICPoint_to_ICHamPoint( parameters, p, ics );

		if ( status.MPI_TAG == EXIT_TAG )
		{
			cout << "Received exit tag." << endl;
			break;
		}

		trajectory.y0[0] = ics.R;
		trajectory.y0[1] = - ics.pR;
		trajectory.y0[2] = ics.theta;
		trajectory.y0[3] = ics.pT;

		// #####################################################
		// p.v0, p.b -- in SI
		double b_integrand = fabs( 2 * M_PI * p.v0 * p.b ); 

		double v0_integrand = (MU_SI/TWO_PI/kT) * sqrt(MU_SI/TWO_PI/kT) * exp( -MU_SI*p.v0*p.v0 / (2.0 * kT)) * 4 * MYPI * p.v0 * p.v0; 

		double stat_weight = b_integrand * v0_integrand;
		//cout << "stat_weight: " << stat_weight << endl;
		// #####################################################
		
		clock_t start = clock();

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
				
			specfunc_value_classical = SPECFUNC_POWERS_OF_TEN * specfunc_coeff * stat_weight * dipfft;
			classical.specfunc_package.push_back( specfunc_value_classical );
	
			spectrum_value_classical = SPECTRUM_POWERS_OF_TEN * spectrum_coeff * stat_weight * omega *  ( 1.0 - exp( - constants::PLANCKCONST_REDUCED * omega / kT ) ) * dipfft;
			classical.spectrum_package.push_back( spectrum_value_classical );

			classical.m2_package += spectrum_value_classical * FREQ_STEP; 
			// PLANCKCONST/2/PI = PLANCKCONST_REDUCED

			d1.correct( classical, omega, FREQ_STEP, kT );
			d2.correct( classical, omega, FREQ_STEP, kT, true );
			d3.correct( classical, omega, FREQ_STEP, kT );
		}

		cout << "(" << world_rank << ") Processing " << p.counter << " trajectory. npoints = " << npoints << "; time = " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "s" << endl;
		
		//// #################################################
		// Sending data
		classical.send();
		d1.send();
		d2.send();
		d3.send();

		classical.clear_package();
		d1.clear_package();
		d2.clear_package();
		d3.clear_package();
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
