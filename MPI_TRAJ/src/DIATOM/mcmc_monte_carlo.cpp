#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>

// checking if dir exists
#include <dirent.h>
#include <sys/stat.h> // mkdir

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

// physical constants
#include "constants.h"

// Gear header files
#include "basis.h"
#include "vmblock.h"
#include "gear.h"
#include "t_dgls.h"

// library for FFT
#include <fftw3.h>

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

void add_vectors( vector<double>& v1, vector<double>& v2 )
{
	assert( v1.size() == v2.size() );

	for ( int i = 0; i < v1.size(); i++ )
	{
		v1[i] += v2[i];
	}
}

void zero_out_vector( vector<double>& v )
{
	for ( int i = 0; i < v.size(); i++ )
	{
		v[i] = 0.0;
	}
}

void mult_vector( vector<double>& v, const double& coeff )
{
	for ( int i = 0; i < v.size(); i++ )
	{
		v[i] *= coeff;
	}
}	

string modify_filename( string filename, string modifier )
{
	size_t dot_position = filename.find( "." );
	if ( dot_position != string::npos )
	{
		string str = filename.substr( 0, dot_position );
		//cout << "str: " << str << endl;
			
		str = str + "_" + modifier + ".txt";

		return str;
	}
	else
	{
		cout << "There is not dot in filename! Setting filename to base" << endl;
		
		return modifier + ".txt";
	}
}

bool check_dir_exists( string dirname )
{
	DIR* dir = opendir( dirname.c_str() );
	if ( dir )
	{
		// directory exists
		closedir( dir );
		return true;
	}
	else if ( ENOENT == errno )
	{
		// directory doesn't exist
		return false;
	}
	else
	{
		// unknown wtf
		cout << "Unknown error status while checking ouput dir" << endl;
		return false;
	}
}

void create_dir( string dirname )
{
	int status = mkdir ( dirname.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
}

void saving_procedure( Parameters& parameters, vector<double>& freqs, vector<double>& total_specfunc, vector<double>& total_spectrum, const double& uniform_integrated_spectrum_total, string modifier = "class" )
{
	string out_dir = parameters.output_directory;
	bool status = check_dir_exists( out_dir );
	if ( status == true )
	{
		//cout << "Directory exists" << endl;
	}
	else
	{
		//cout << "Creating directory" << endl;
		create_dir( out_dir );
	}

	if ( modifier == "class" )
	{
		save( freqs, total_specfunc, out_dir + "/" + parameters.specfunc_filename );
		save( freqs, total_spectrum, out_dir + "/" + parameters.spectrum_filename );
		save( uniform_integrated_spectrum_total, out_dir + "/" + parameters.m2_filename );
	}
	else	
	{
		string specfunc_filename = modify_filename( parameters.specfunc_filename, modifier );
		save( freqs, total_specfunc, out_dir + "/" + specfunc_filename );

		string spectrum_filename = modify_filename( parameters.spectrum_filename, modifier );
		save( freqs, total_spectrum, out_dir + "/" + spectrum_filename );

		string m2_filename = modify_filename( parameters.m2_filename, modifier );
		save( uniform_integrated_spectrum_total, out_dir + "/" + m2_filename );
	}
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
	
	// ##########################################################
	vector<double> specfunc_package_class( FREQ_SIZE );
	vector<double> spectrum_package_class( FREQ_SIZE );
	vector<double> chunk_specfunc_class( FREQ_SIZE );
	vector<double> chunk_spectrum_class( FREQ_SIZE );
	vector<double> total_spectrum_class( FREQ_SIZE );
	vector<double> total_specfunc_class( FREQ_SIZE );

	double uniform_integrated_spectrum_package_class;
	double uniform_integrated_spectrum_chunk_class = 0.0;
	double uniform_integrated_spectrum_total_class = 0.0;
	// ##########################################################

	// ##########################################################
	vector<double> spectrum_package_d1( FREQ_SIZE );
	vector<double> specfunc_package_d1( FREQ_SIZE );
	vector<double> chunk_spectrum_d1( FREQ_SIZE );
	vector<double> chunk_specfunc_d1( FREQ_SIZE );
	vector<double> total_spectrum_d1( FREQ_SIZE );
	vector<double> total_specfunc_d1( FREQ_SIZE );

	double uniform_integrated_spectrum_package_d1;
	double uniform_integrated_spectrum_chunk_d1 = 0.0;
	double uniform_integrated_spectrum_total_d1 = 0.0;
	// ##########################################################

	// ##########################################################
	vector<double> spectrum_package_d2( FREQ_SIZE );
	vector<double> specfunc_package_d2( FREQ_SIZE );
	vector<double> chunk_spectrum_d2( FREQ_SIZE );
	vector<double> chunk_specfunc_d2( FREQ_SIZE );
	vector<double> total_spectrum_d2( FREQ_SIZE );
	vector<double> total_specfunc_d2( FREQ_SIZE );

	double uniform_integrated_spectrum_package_d2;
	double uniform_integrated_spectrum_chunk_d2 = 0.0;
	double uniform_integrated_spectrum_total_d2 = 0.0;
	// ##########################################################

	// ##########################################################
	vector<double> spectrum_package_d3( FREQ_SIZE );
	vector<double> specfunc_package_d3( FREQ_SIZE );
	vector<double> chunk_spectrum_d3( FREQ_SIZE );
	vector<double> chunk_specfunc_d3( FREQ_SIZE );
	vector<double> total_spectrum_d3( FREQ_SIZE );
	vector<double> total_specfunc_d3( FREQ_SIZE );

	double uniform_integrated_spectrum_package_d3;
	double uniform_integrated_spectrum_chunk_d3 = 0.0;
	double uniform_integrated_spectrum_total_d3 = 0.0;
	// ##########################################################

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

	//// sending first trajectory	
	//for ( int i = 1; i < world_size; i++ )
	//{
		//p = parameters.generate_uniform_point( B_CHUNKS_VECTOR[ b_chunk_counter ],
												   //B_CHUNKS_VECTOR[ b_chunk_counter + 1 ],
												   //V0_CHUNKS_VECTOR[ v0_chunk_counter ],
												   //V0_CHUNKS_VECTOR[ v0_chunk_counter + 1 ]);
		//p.counter = sent; 
			
		////cout << "###" << endl;
		////cout << "generated p.b: " << p.b << endl;
		////cout << "generated p.v0: " << p.v0 << endl;
		////cout << "p.counter: " << p.counter << endl;
		////cout << "###" << endl;

		//MPI_Send( &p, 1, MPI_ICPoint, i, 0, MPI_COMM_WORLD );
		//sent++;
	//}

	//clock_t start = clock();

	//while( true )
	//{
		//if ( is_finished )	
		//{
			//for ( int i = 1; i < world_size; i++ )
			//{	
				//MPI_Send( &is_finished, 1, MPI_INT, i, EXIT_TAG, MPI_COMM_WORLD );
			//}

			//break;
		//}
		
		//// ############################################################
		//// Receiving data
		//MPI_Recv( &specfunc_package_class[0], FREQ_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		//source = status.MPI_SOURCE;

		//MPI_Recv( &spectrum_package_class[0], FREQ_SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		//MPI_Recv (&uniform_integrated_spectrum_package_class, 1, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

		//received++;
		
		//if ( parameters.d1_status == true )
		//{
			//MPI_Recv( &specfunc_package_d1[0], FREQ_SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//MPI_Recv( &spectrum_package_d1[0], FREQ_SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//MPI_Recv( &uniform_integrated_spectrum_package_d1, 1, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		//}

		//if ( parameters.d2_status == true )
		//{
			//MPI_Recv( &specfunc_package_d2[0], FREQ_SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//MPI_Recv( &spectrum_package_d2[0], FREQ_SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//MPI_Recv( &uniform_integrated_spectrum_package_d2, 1, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		//}

		//if ( parameters.d3_status == true )
		//{
			//MPI_Recv( &specfunc_package_d3[0], FREQ_SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//MPI_Recv( &spectrum_package_d3[0], FREQ_SIZE, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			//MPI_Recv( &uniform_integrated_spectrum_package_d3, 1, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		//}

		//// ############################################################

		////cout << "sent: " << sent << "; received: " << received << "; NPOINTS: " << parameters.NPOINTS << endl;	
		//if ( received == parameters.NPOINTS )
		//{
			//// multiplyting spectrum and spectral function by AREA and dividing by NPOINTS
			//B_CHUNK = B_CHUNKS_VECTOR[b_chunk_counter + 1] -
					  //B_CHUNKS_VECTOR[b_chunk_counter];
			//V0_CHUNK = V0_CHUNKS_VECTOR[v0_chunk_counter + 1] -
					   //V0_CHUNKS_VECTOR[v0_chunk_counter];
			
			//double multiplier = B_CHUNK * V0_CHUNK / parameters.NPOINTS; 
			////cout << "multiplier: " << multiplier << endl;

			//// ##################################################################################
			//mult_vector( chunk_specfunc_class, multiplier );
			//mult_vector( chunk_spectrum_class, multiplier );
			//uniform_integrated_spectrum_chunk_class *= multiplier; 

			//if ( parameters.d1_status == true )
			//{
				//mult_vector( chunk_specfunc_d1, multiplier );
				//mult_vector( chunk_spectrum_d1, multiplier );

				//uniform_integrated_spectrum_chunk_d1 *= multiplier;
			//}
			
			//if ( parameters.d2_status == true )
			//{
				//mult_vector( chunk_specfunc_d2, multiplier );
				//mult_vector( chunk_spectrum_d2, multiplier );

				//uniform_integrated_spectrum_chunk_d2 *= multiplier;
			//}

			//if ( parameters.d3_status == true )
			//{
				//mult_vector( chunk_specfunc_d3, multiplier );
				//mult_vector( chunk_spectrum_d3, multiplier );

				//uniform_integrated_spectrum_chunk_d3 *= multiplier;
			//}

			//// ##################################################################################

			//cout << "#####################################" << endl;
			//cout << "Current chunk: (b) " << b_chunk_counter << " (v0) " << v0_chunk_counter << endl;
			//cout << "B: " << B_CHUNKS_VECTOR[b_chunk_counter] << " -- " << B_CHUNKS_VECTOR[b_chunk_counter + 1] << endl;
			//cout << "V0: " << V0_CHUNKS_VECTOR[v0_chunk_counter] << " -- " << V0_CHUNKS_VECTOR[v0_chunk_counter + 1] << endl;
			//cout << "M2 (class): " << uniform_integrated_spectrum_chunk_class << endl;
			//cout << "M2 (d1): " << uniform_integrated_spectrum_chunk_d1 << endl;
			//cout << "M2 (d2): " << uniform_integrated_spectrum_chunk_d2 << endl;
			//cout << "M2 (d3): " << uniform_integrated_spectrum_chunk_d3 << endl;
			//cout << "Time for chunk: " << (clock() - start) / (double) CLOCKS_PER_SEC << "s" << endl;
			//start = clock();
			//cout << "#####################################" << endl;
			
			//// adding chunk spectrum to total; zeroing out chunk data
			//cout << ">> Added chunk spectrum to total" << endl << endl;	
			//add_vectors( total_specfunc_class, chunk_specfunc_class );
			//add_vectors( total_spectrum_class, chunk_spectrum_class );
			
			//zero_out_vector( chunk_specfunc_class );
			//zero_out_vector( chunk_spectrum_class );
			
			//uniform_integrated_spectrum_total_class += uniform_integrated_spectrum_chunk_class;
			//uniform_integrated_spectrum_chunk_class = 0.0;
			
			//if ( parameters.d1_status == true )
			//{
				//add_vectors( total_specfunc_d1, chunk_specfunc_d1 );
				//add_vectors( total_spectrum_d1, chunk_spectrum_d1 );

				//zero_out_vector( chunk_specfunc_d1 );
				//zero_out_vector( chunk_spectrum_d1 );

				//uniform_integrated_spectrum_total_d1 += uniform_integrated_spectrum_chunk_d1;
				//uniform_integrated_spectrum_chunk_d1 = 0.0;
			//}
		
			//if ( parameters.d2_status == true )
			//{
				//add_vectors( total_specfunc_d2, chunk_specfunc_d2 );
				//add_vectors( total_spectrum_d2, chunk_spectrum_d2 );

				//zero_out_vector( chunk_specfunc_d2 );
				//zero_out_vector( chunk_spectrum_d2 );

				//uniform_integrated_spectrum_total_d2 += uniform_integrated_spectrum_chunk_d2;
				//uniform_integrated_spectrum_chunk_d2 = 0.0;
			//}
			
			//if ( parameters.d3_status == true )
			//{
				//add_vectors( total_specfunc_d3, chunk_specfunc_d3 );
				//add_vectors( total_spectrum_d3, chunk_spectrum_d3 );

				//zero_out_vector( chunk_specfunc_d3 );
				//zero_out_vector( chunk_spectrum_d3 );

				//uniform_integrated_spectrum_total_d3 += uniform_integrated_spectrum_chunk_d3;
				//uniform_integrated_spectrum_chunk_d3 = 0.0;
			//}

			////cout << "point counter is set to 0" << endl << endl;
			//sent = 0;
			//received = 0;
			//// ##################################################

			//cout << "Saving spectrum" << endl;
			//saving_procedure( parameters, freqs, total_specfunc_class, total_spectrum_class, uniform_integrated_spectrum_total_class ); 
	
			//if ( parameters.d1_status == true )
			//{
				//saving_procedure( parameters, freqs, total_specfunc_d1, total_spectrum_d1, uniform_integrated_spectrum_total_d1, "d1" );
			//}

			//if ( parameters.d2_status == true )
			//{
				//saving_procedure( parameters, freqs, total_specfunc_d2, total_spectrum_d2, uniform_integrated_spectrum_total_d2, "d2" );
			//}

			//if ( parameters.d3_status == true )
			//{
				//saving_procedure( parameters, freqs, total_specfunc_d3, total_spectrum_d3, uniform_integrated_spectrum_total_d3, "d3" );
			//}

			//// ##################################################
			//if ( b_chunk_counter < b_chunk_max - 1 )
			//{
				//cout << ">> b_chunk_counter is not maximum. incrementing..." << endl;
				//b_chunk_counter++;
			//}	
			//else
			//{
				//if ( v0_chunk_counter < v0_chunk_max - 1 )
				//{
					//cout << ">> v0_chunk_counter is not maximum. Incrementing..." << endl;
					//cout << ">> setting b_chunk_counter to 0" << endl;
					//v0_chunk_counter++;
					//b_chunk_counter = 0;
				//}
				//else
				//{
					//cout << ">> all chunks are integrated!" << endl;
					//// setting finished flag to true
					//is_finished = true;
					//// starting new iteration which checks the finished flag 
					//continue;
				//}
			//}
			//// ##################################################

			//// ##################################################
			//for ( int i = 1; i < world_size; i++ )
			//{
				//p = parameters.generate_uniform_point( B_CHUNKS_VECTOR[ b_chunk_counter ],
												   	   //B_CHUNKS_VECTOR[ b_chunk_counter + 1 ],
												   	   //V0_CHUNKS_VECTOR[ v0_chunk_counter ],
													   //V0_CHUNKS_VECTOR[ v0_chunk_counter + 1 ]);
				//p.counter = sent; 	
			
				////cout << "###" << endl;
				////cout << "generated p.b: " << p.b << endl;
				////cout << "generated p.v0: " << p.v0 << endl;
				////cout << "p.counter: " << p.counter << endl;
				////cout << "###" << endl;

				//MPI_Send( &p, 1, MPI_ICPoint, i, 0, MPI_COMM_WORLD );
				////cout << "point_counter: " << point_counter << endl;
				//sent++;
			//}
			//// ##################################################
		//}		
		
		//// ############################################################
		//for ( int i = 0; i < FREQ_SIZE; i++ )
		//{
			//chunk_specfunc_class[i] += specfunc_package_class[i];
			//chunk_spectrum_class[i] += spectrum_package_class[i];
		//}
		//uniform_integrated_spectrum_chunk_class += uniform_integrated_spectrum_package_class;

		//if ( parameters.d1_status == true )
		//{
			//for ( int i = 0; i < FREQ_SIZE; i++ )
			//{
				//chunk_specfunc_d1[i] += specfunc_package_d1[i];
				//chunk_spectrum_d1[i] += spectrum_package_d1[i];
			//}
			//uniform_integrated_spectrum_chunk_d1 += uniform_integrated_spectrum_package_d1;
		//}
		
		//if ( parameters.d2_status == true )
		//{
			//for ( int i = 0; i < FREQ_SIZE; i++ )
			//{
				//chunk_specfunc_d2[i] += specfunc_package_d2[i];
				//chunk_spectrum_d2[i] += spectrum_package_d2[i];
			//}
			//uniform_integrated_spectrum_chunk_d2 += uniform_integrated_spectrum_package_d2;
		//}
		
		//if ( parameters.d3_status == true )
		//{
			//for ( int i = 0; i < FREQ_SIZE; i++ )
			//{
				//chunk_specfunc_d3[i] += specfunc_package_d3[i];
				//chunk_spectrum_d3[i] += spectrum_package_d3[i];
			//}
			//uniform_integrated_spectrum_chunk_d3 += uniform_integrated_spectrum_package_d3;
		//}
		
		////cout << "Added package spectrum to chunk." << endl;
		////cout << "chunk_specfunc[0]: " << chunk_specfunc[0] << endl;
		//// ############################################################
		
		
		//if ( sent < parameters.NPOINTS )
		//{
			//p = parameters.generate_uniform_point( B_CHUNKS_VECTOR[ b_chunk_counter ],
												   //B_CHUNKS_VECTOR[ b_chunk_counter + 1 ],
												   //V0_CHUNKS_VECTOR[ v0_chunk_counter ],
												   //V0_CHUNKS_VECTOR[ v0_chunk_counter + 1 ]);
			//p.counter = sent; 
			
			////cout << "###" << endl;
			////cout << "generated p.b: " << p.b << endl;
			////cout << "generated p.v0: " << p.v0 << endl;
			////cout << "p.counter: " << p.counter << endl;
			////cout << "###" << endl;

			//MPI_Send( &p, 1, MPI_ICPoint, source, 0, MPI_COMM_WORLD );
			////cout << "point_counter: " << point_counter << endl;
			//sent++;
		//}
	//}
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

	//while ( true )
	//{
		//MPI_Recv( &p, 1, MPI_ICPoint, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status ); 
		////cout << "(" << world_rank << ") Received ICPoint: " << p.counter << endl;
		////show_point( p );
		//transform_ICPoint_to_ICHamPoint( parameters, p, ics );

		//if ( status.MPI_TAG == EXIT_TAG )
		//{
			//cout << "Received exit tag." << endl;
			//break;
		//}

		//y0[0] = ics.R;
		//y0[1] = - ics.pR;
		//y0[2] = ics.theta;
		//y0[3] = ics.pT;

		////cout << "####" << endl;
		////cout << "p.v0: " << p.v0 << endl;
		////cout << "p.b: " << p.b << endl;
		////cout << "ics.R = " << y0[0] << endl;
		////cout << "ics.pR = " << y0[1] << endl;
		////cout << "ics.theta = " << y0[2] << endl;
		////cout << "ics.pT = " << y0[3] << endl;
		////cout << "#####" << endl;

		//// out of memory?
		//if ( !vmcomplete(vmblock) )
		//{ 
			//cout << "mgear: out of memory" << endl;
			//return;
		//}

		//// #####################################################
		//// p.v0, p.b -- in SI
		//double b_integrand = fabs( 2 * M_PI * p.v0 * p.b ); 

		//double v0_integrand = (MU_SI/TWO_PI/kT) * sqrt(MU_SI/TWO_PI/kT) * exp( -MU_SI*p.v0*p.v0 / (2.0 * kT)) * 4 * MYPI * p.v0 * p.v0; 

		//double stat_weight = b_integrand * v0_integrand;
		////cout << "stat_weight: " << stat_weight << endl;
		//// #####################################################
		
		//int counter = 0;
		//double R_end_value = y0[0] + 0.001;

		//// dipole moment in laboratory frame
		//vector<double> temp( 3 );
		//vector<double> dipx;
		//vector<double> dipy;
		//vector<double> dipz;
		
		//// #####################################################
		//t0 = 0.0;

		//h = 0.1;    		  // initial step size
		//xend = parameters.sampling_time; // initial right bound of integration
		//// #####################################################

		//clock_t start = clock();

		//// #####################################################
		//while ( y0[0] < R_end_value ) 
		//{
			//if ( counter == parameters.MaxTrajectoryLength )
			//{
				////cout << "Trajectory cut!" << endl;
				//break;
			//}
			
			//fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);
			
			////cout << "%%%" << endl;
			////cout << "t0: " << t0 << endl;
			////cout << "xend: " << xend << endl;
			////cout << "%%%" << endl;

			//if ( fehler != 0 ) 
			//{
				//cout << "Gear4: error n = " << 10 + fehler << endl;
				//break;
			//}

			////cout << "t0: " << t0 << "; r: " << y0[0] << endl;

			//// y0[0] -- R
			//// y0[1] -- PR
			//// y0[2] -- \theta
			//// y0[3] -- p_\theta
			
			//if ( parameters.use_S_matrix == true )
			//{
				//transform_dipole( temp, y0[0], y0[2] );
			//}
			//else
			//{
				//dipole_without_S( temp, y0[0] );
				////cout << "transformed dipole" << endl;
			//}
			
			////cout << "t: " << t0 * constants::ATU <<
				   	////"; R (alu): " << y0[0] << 	
					////"; R: " << y0[0] * constants::ALU << 
					////"; dipole z: " << temp[2] * constants::ADIPMOMU << endl;

			//dipx.push_back( temp[0] );
			//dipy.push_back( temp[1] );
			//dipz.push_back( temp[2] );

			//xend = parameters.sampling_time * (counter + 2);

			//aufrufe = 0;  // actual number of calls

			//counter++;
		//}
		//// #####################################################
	
		//// #####################################################
		//// length of dipole vector = number of samples
		//int npoints = dipz.size();

		////cout << "FREQ_SIZE: " << FREQ_SIZE << endl;

		//// zeroing input arrays
		//fourier.zero_out_input( );

		//// copying data to them
		//copy_to( dipx, fourier.inx );
		//copy_to( dipy, fourier.iny );
		//copy_to( dipz, fourier.inz );

		//// executing fourier transform
		//fourier.do_fourier( );

		//double omega, dipfft;

		//double ReFx, ReFy, ReFz;
		//double ImFx, ImFy, ImFz;
		
		//vector<double> specfunc_class;
		//vector<double> spectrum_class;

		//vector<double> specfunc_d1;
		//vector<double> specfunc_d2;
		//vector<double> specfunc_d3;
		//vector<double> specfunc_d4;

		//vector<double> spectrum_d1;
		//vector<double> spectrum_d2;
		//vector<double> spectrum_d3;
		//vector<double> spectrum_d4;

		//double specfunc_value_class;
		//double spectrum_value_class;

		//double uniform_integrated_spectrum_class = 0;
		//double uniform_integrated_spectrum_d1 = 0;
		//double uniform_integrated_spectrum_d2 = 0;
		//double uniform_integrated_spectrum_d3 = 0;
		//double uniform_integrated_spectrum_d4 = 0;

		//for ( int k = 0; k < FREQ_SIZE; k++ )
		//{
				//omega = 2.0 * M_PI * constants::LIGHTSPEED_CM * freqs[k];

				//ReFx = fourier.outx[k][0];
				//ReFy = fourier.outy[k][0];
				//ReFz = fourier.outz[k][0];
				//ImFx = fourier.outx[k][1];
				//ImFy = fourier.outy[k][1];
				//ImFz = fourier.outz[k][1];

				//dipfft = ReFx * ReFx + ReFy * ReFy + ReFz * ReFz +
						 //ImFx * ImFx + ImFy * ImFy + ImFz * ImFz;
				////cout << "dipfft[" << k << "] = " << dipfft * constants::ADIPMOMU * constants::ADIPMOMU << endl;	

				//specfunc_value_class = SPECFUNC_POWERS_OF_TEN * specfunc_coeff * stat_weight * dipfft;
				//specfunc_class.push_back( specfunc_value_class );	

				//spectrum_value_class = SPECTRUM_POWERS_OF_TEN * spectrum_coeff * stat_weight * omega *  ( 1.0 - exp( - constants::PLANCKCONST_REDUCED * omega / kT ) ) * dipfft;
				//spectrum_class.push_back( spectrum_value_class );

				//uniform_integrated_spectrum_class += spectrum_value_class * FREQ_STEP; 
				//// PLANCKCONST/2/PI = PLANCKCONST_REDUCED

				//if ( parameters.d1_status == true )
				//{
					//double d1 = 2.0 / (1.0 + exp(-constants::PLANCKCONST_REDUCED * omega / kT));
					//double specfunc_value_d1 = specfunc_value_class * d1;
					//double spectrum_value_d1 = spectrum_value_class * d1;
										
					//specfunc_d1.push_back( specfunc_value_d1 );
					//spectrum_d1.push_back( spectrum_value_d1 );

					//uniform_integrated_spectrum_d1 += spectrum_value_d1 * FREQ_STEP;
				//}

				//if ( parameters.d2_status == true )
				//{
					//double d2;
					//if ( omega == 0 ) d2 = 0.0;	
					//else d2 = constants::PLANCKCONST_REDUCED * omega / kT / (1.0 - exp(-constants::PLANCKCONST_REDUCED * omega / kT));

					//double specfunc_value_d2 = specfunc_value_class * d2;
					//double spectrum_value_d2 = spectrum_value_class * d2;
					
					//specfunc_d2.push_back( specfunc_value_d2 );
					//spectrum_d2.push_back( spectrum_value_d2 );
							
					//uniform_integrated_spectrum_d2 += spectrum_value_d2 * FREQ_STEP;
				//}

				//if ( parameters.d3_status == true )
				//{
					//double d3 = exp( constants::PLANCKCONST_REDUCED * omega / kT / 2.0 );
				   	//double specfunc_value_d3 = specfunc_value_class * d3;
					//double spectrum_value_d3 = spectrum_value_class * d3;

					//specfunc_d3.push_back( specfunc_value_d3 );
					//spectrum_d3.push_back( spectrum_value_d3 );

					//uniform_integrated_spectrum_d3 += spectrum_value_d3 * FREQ_STEP;
				//}

				//if ( parameters.d4_status == true )
				//{
					//// to be done
				//}
		//}

		//cout << "(" << world_rank << ") Processing " << p.counter << " trajectory. npoints = " << npoints << "; time = " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "s" << endl;
		
		////// #################################################
		//// Sending data
		//MPI_Send( &specfunc_class[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		//MPI_Send( &spectrum_class[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		//MPI_Send( &uniform_integrated_spectrum_class, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );	
		
		//if ( parameters.d1_status == true )
		//{
			//MPI_Send( &specfunc_d1[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			//MPI_Send( &spectrum_d1[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			//MPI_Send( &uniform_integrated_spectrum_d1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		//}

		//if ( parameters.d2_status == true )
		//{
			//MPI_Send( &specfunc_d2[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			//MPI_Send( &spectrum_d2[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			//MPI_Send( &uniform_integrated_spectrum_d2, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		//}

		//if ( parameters.d3_status == true )
		//{
			//MPI_Send( &specfunc_d3[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			//MPI_Send( &spectrum_d3[0], FREQ_SIZE, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			//MPI_Send( &uniform_integrated_spectrum_d3, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		//}
		//// #################################################
	//}
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
