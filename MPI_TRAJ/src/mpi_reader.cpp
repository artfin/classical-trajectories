#include <mpi.h>

#include <iostream>
#include <fstream>

// string conversion
#include <string>
#include <sstream>

#include <vector>
#include <chrono>

using namespace std;

// exiting tag
const int EXIT_TAG = 42;

// output directory to process
//const string OUTPUT_DIR = "output/dips/";
const string OUTPUT_DIR = "first_exp/dips/";

// frequency bins bounds
const double MIN_FREQ = 0.0;
const double MAX_FREQ = 200.0;
const double BIN_SIZE = 0.5;

// number of bins
const int NBINS = ( MAX_FREQ - MIN_FREQ ) / BIN_SIZE;

int main( int argc, char* argv[] )
{
	MPI_Init( &argc, &argv );

	int world_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );

	int world_size;
	MPI_Comm_size( MPI_COMM_WORLD, &world_size );

	// variable to store status of messages
	MPI_Status status;
	int source;

	// root 	
	if ( world_rank == 0 )
	{
		cout << "argc: " << argc << endl;
		
		for ( int i = 0; i < argc; i++ )
		{
			cout << "argv[" << i << "]: " << argv[i] << endl;
		}

		cout << "---------------------------------------------" << endl;

		if ( argc != 3 )
		{
			cerr << "USAGE: ./... (int) starting file number (int) ending file number" << endl;
		   	MPI_Finalize();
			
			exit( 1 );	
		}

		int START_filenumber = atoi( argv[1] );
		int END_filenumber = atoi( argv[2] );
		
		int curr_filenumber = START_filenumber;

		// right bounds of frequency bins
		vector<double> freqs_rbound( NBINS );
		vector<double> intensities( NBINS );

		// prefilling frequency vector
		for ( int i = 0; i < freqs_rbound.size(); i++ )
		{
			freqs_rbound[i] = BIN_SIZE * i; 
		}

		// tracking time
		auto startTime = chrono::high_resolution_clock::now();

		// sending first message to slaves
		for ( int i = 1; i < world_size; i++ )
		{
			//cout << "Root sends curr_filenumber = " << curr_filenumber << " to process " << i << endl;

			MPI_Send( &curr_filenumber, 1, MPI_INT, i, 0, MPI_COMM_WORLD );
		   	curr_filenumber++;	
		}

		// number of slaves alive 
		int alive = world_size - 1;
		
		while ( true )
		{
			if ( curr_filenumber % 1000 == 0 )
			{
				auto partial = chrono::high_resolution_clock::now();
				cout << "Curr_filenumber: " << curr_filenumber << 
						"; time elapsed: " << chrono::duration_cast<chrono::milliseconds> ( partial - startTime ).count() / 1000.0 << "s" << endl;
			}

			if ( alive == 0 )
			{
				break;
			}
	
			// size of vector to get
			int package_size; 
			MPI_Recv( &package_size, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		   	
			source = status.MPI_SOURCE;
			//cout << ">> Root received report from process " << source << endl;
			
			vector<double> freqs_package( package_size );
			vector<double> intensities_package( package_size );
			
			// file is not empty
			if ( package_size != 0 )
			{
				MPI_Recv( &freqs_package[0], package_size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
				//cout << ">> Root received frequency package of size " << package_size << endl;

				MPI_Recv( &intensities_package[0], package_size, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
				//cout << ">> Root received intensities package of size " << package_size << endl;
			}

			for ( int i = 0; i < package_size; i++ )
			{
				for ( int j = 0; j < NBINS - 1; j++ )
				{
					if ( freqs_package[i] >= freqs_rbound[j] && freqs_package[i] < freqs_rbound[j + 1] )
					{
						intensities[j] += intensities_package[i];
					}	
				}
			}

			//cout << ">> Root've done binning." << endl;

			// if there is file to be processed
			if ( curr_filenumber <= END_filenumber )
			{
				MPI_Send( &curr_filenumber, 1, MPI_INT, source, 0, MPI_COMM_WORLD );
				//cout << "Root send curr_filenumber = " << curr_filenumber << " to process " << source << endl;
				curr_filenumber++;
			}
			// if not -- killing slave
			else
			{
				MPI_Send( &curr_filenumber, 1, MPI_INT, source, EXIT_TAG, MPI_COMM_WORLD );
				//cout << ">> Root sends killing message to process " << source << endl; 
				// decrementing alive processes counter
				alive--;
			}
		}

		auto endTime = chrono::high_resolution_clock::now();
		cout << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds> ( endTime - startTime ).count() / 1000.0 << "s" << endl;;


		ofstream output_file;
		output_file.open( "spectra.dat" );

		for ( int i = 0; i < freqs_rbound.size(); i++ )
		{
			output_file << freqs_rbound[i] << "  " << intensities[i] << endl;
		}

		output_file.close();
	}
	// slave
	else
	{
		// current file to process
		int curr_filenumber;

		while ( true )
		{
			// receiving number of file to process
			MPI_Recv( &curr_filenumber, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if ( status.MPI_TAG == EXIT_TAG )
			{
				//cout << "Slave process " << world_rank << " receives kill message." << endl;
				break;
			}	

			//cout << "Slave process: " << world_rank << "; received curr_filenumber = " << curr_filenumber << " from root" << endl;

			ostringstream strs;
			strs << curr_filenumber;
			string filename = OUTPUT_DIR + strs.str() + ".bin";

			ifstream file( filename, ios::in | ios::binary | ios::ate );
			
			streampos file_size = file.tellg();
			int ndoubles = file_size / 8;

			char* memblock = new char [file_size];
			file.seekg(0, ios::beg );
			file.read( memblock, file_size );
			file.close();

			double* arr = (double*) memblock;
		
			int package_size;	
			
			// if there is nothing to process
			if ( ndoubles == 0 )
			{
				package_size = 0;
				MPI_Send( &package_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
			}
			else
			{
				vector<double> freqs( ndoubles / 2 );
				vector<double> intensities( ndoubles / 2 );
				
				for( int cycle_counter = 0, vec_counter = 0; cycle_counter < ndoubles; cycle_counter += 2, vec_counter++ )
				{
					freqs[vec_counter] = arr[cycle_counter];
					intensities[vec_counter] = arr[cycle_counter + 1];
				}

				//for ( int i = 0; i < freqs.size(); i++ )
				//{
						//cout << "ON SLAVE i: " << i << 
								//"; freqs: " << freqs[i] << 
								//"; intensities: " << intensities[i] << endl;
				//}

				// sending read data
				package_size = ndoubles / 2; 
				MPI_Send( &package_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
			
				//cout << "Slave sending packages of size " << package_size << endl; 
				MPI_Send( &freqs[0], package_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
				MPI_Send( &intensities[0], package_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );

			}
		}
	}

	MPI_Finalize();

	return 0;
}
