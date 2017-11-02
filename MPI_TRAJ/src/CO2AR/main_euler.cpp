#include <mpi.h>

#include <iostream>

// matrix multiplication
#include "matrix_euler.h"

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

// gsl histogram
#include <gsl/gsl_histogram.h>

// ############################################
// Exit tag for killing slave
const int EXIT_TAG = 42;
// ############################################

// ############################################
// Number of initial conditions per trajectory
const int ICPERTRAJ = 11; 
// ############################################

// ############################################
const double Temperature = 298; 
// ############################################

// ############################################
// GSL Histogram parameters (spectrum boundaries)
const int NBINS = 400;
const double LBOUND = 0.0;
const double RBOUND = 400.0;
// ############################################

using namespace std;

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
  	
  	double *out = new double[10];

	//cout << "###########" << endl;
	//cout << "inside syst" << endl;
	//cout << "y[0]: " << y[0] << endl;
	//cout << "y[1]: " << y[1] << endl;
	//cout << "y[2]: " << y[2] << endl;
	//cout << "y[3]: " << y[3] << endl;
	//cout << "y[4]: " << y[4] << endl;
	//cout << "y[5]: " << y[5] << endl;
	//cout << "y[6]: " << y[6] << endl;
	//cout << "y[7]: " << y[7] << endl;
	//cout << "y[8]: " << y[8] << endl;
	//cout << "y[9]: " << y[9] << endl;
	//cout << "#############" << endl;

  	rhs(out, y[0], y[1], // R Theta 
			 y[2], y[3], // pR pTheta
		   	 y[4], y[5], y[6], // phi theta psi
		     y[7], y[8], y[9] // p_phi p_theta p_psi
	    );

	//cout << "out[0] = " << out[0] << endl;
	//cout << "out[1] = " << out[1] << endl;
	//cout << "out[2] = " << out[2] << endl;
	//cout << "out[3] = " << out[3] << endl;
	//cout << "out[4] = " << out[4] << endl;
	//cout << "out[5] = " << out[5] << endl;
	//cout << "out[6] = " << out[6] << endl;
	//cout << "out[7] = " << out[7] << endl;
	//cout << "out[8] = " << out[8] << endl;
	//cout << "out[9] = " << out[9] << endl;

    f[0] = out[0]; // dR/dt  
    f[1] = out[1]; // d(Theta)/dt
    f[2] = out[2]; // d(pR)/dt
    f[3] = out[3]; // d(pT)/dt
    f[4] = out[4]; // d(phi)/dt 
    f[5] = out[5]; // d(theta)/dt
    f[6] = out[6]; // d(psi)/dt
    f[7] = out[7]; // d(p_phi)/dt
    f[8] = out[8]; // d(p_theta)/dt
    f[9] = out[9]; // d(p_psi)/dt

    delete [] out;
}

void master_code( int world_size )
{
	MPI_Status status;
	int source;

	FILE* inputfile = fopen("input/CO2AR/ics_euler", "r");

	// counter of calculated trajectories
	int NTRAJ = 0;

	double *ics = new double[ICPERTRAJ];

	// result of reading values from file
	int scanfResult;

	// sending first message to slaves
	for ( int i = 1; i < world_size; i++ ) 
	{
		scanfResult = fwscanf(inputfile, L"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &ics[0], &ics[1], &ics[2], &ics[3], &ics[4], &ics[5], &ics[6], &ics[7], &ics[8], &ics[9], &ics[10]);
	
		cout << "Read initial condition: " << ics[0] << " " 
										   << ics[1] << " "
										   << ics[2] << " "
										   << ics[3] << " "
										   << ics[4] << " "
										   << ics[5] << " "
										   << ics[6] << " "
										   << ics[7] << " "
										   << ics[8] << " "
										   << ics[9] << " "
										   << ics[10] << endl;
		
		MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);

		NTRAJ++;
	}
	
	// number of alive processes
	int alive = world_size - 1;

	// prepare bins
	gsl_histogram *histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( histogram, LBOUND, RBOUND );

	while ( true )
	{	
		// exit when all slaves are killed
		if ( alive == 0 )
		{
			break;
		}
	
		if ( NTRAJ % 1000 == 0 )
		{
			cout << ">> Saving histogram... " << endl;
			//save_histogram( histogram, NBINS, spectrum_filename );
		}

		// receiving message from any of slaves	
		int package_size;
		MPI_Recv( &package_size, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		source = status.MPI_SOURCE;
	
		cout << "Master received a package size = " << package_size << " from process " << source << endl;

		vector<double> freqs_package( package_size );
		vector<double> intensities_package( package_size );

		if ( package_size != 0 )
		{
			MPI_Recv( &freqs_package[0], package_size, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			cout << "Master received frequency package from process " << source << endl;
			
			MPI_Recv( &intensities_package[0], package_size, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			cout << "Master received intensities_package from process " << source << endl;
		}
		
		// applying immediate binning of values
		for ( int i = 0; i < freqs_package.size(); i++ )
		{
			// increases the value of appropriate bin by intensity
			gsl_histogram_accumulate( histogram, freqs_package[i], intensities_package[i] );
		}

		// reading another line from file
		scanfResult = fwscanf(inputfile, L"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &ics[0], &ics[1], &ics[2], &ics[3], &ics[4], &ics[5], &ics[6], &ics[7], &ics[8], &ics[9], &ics[10]);

		// if it's not ended yet then sending new chunk of work to slave
		if ( scanfResult != -1)
   		{
			cout << "Master sends new initial conditions to process " << source << endl;
			MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);
			NTRAJ++;
		}
		// work is done, sending a killing message
		else
		{
			MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, source, EXIT_TAG, MPI_COMM_WORLD);
			alive--;
		}
	}

	gsl_histogram_free( histogram );

	fftw_cleanup();

	delete [] ics;	

	fclose(inputfile);
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
		
	while ( true )
	{
		MPI_Recv(&ics[0], ICPERTRAJ, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		cout << ">> Process " << world_rank << " received new initial conditions." << endl; 

		if ( status.MPI_TAG == EXIT_TAG )
		{
			break;
		}

		N = 10;
		vmblock = vminit();
		y0 = (REAL*) vmalloc(vmblock, VEKTOR, N, 0);
			
		// out of memory?
  		if ( !vmcomplete(vmblock) )
		{ 
   			printf("mgear: out of memory.\n");
   			return;
 		}

		// according to Ivanov:
		// sampling time = 50 fs

		// in atomic time units
		const double sampling_time = 10.0; 

  		epsabs = 1E-13;
  		epsrel = 1E-13;
  
  		t0 = 0.0;
			
  		h = 0.1;         		// initial step size
  		xend = sampling_time;   // initial right bound of integration
  		fmax = 1e8;  	 		// maximal number of calls 
			
		//double *dipole = new double [3];
		copy_initial_conditions( ics, y0, N );

		cout << "y0[0] = " << y0[0] << endl;
		cout << "y0[1] = " << y0[1] << endl;		
		cout << "y0[2] = " << y0[2] << endl;		
		cout << "y0[3] = " << y0[3] << endl;		
		cout << "y0[4] = " << y0[4] << endl;		
		cout << "y0[5] = " << y0[5] << endl;		
		cout << "y0[6] = " << y0[6] << endl;		
		cout << "y0[7] = " << y0[7] << endl;		
		cout << "y0[8] = " << y0[8] << endl;		
		cout << "y0[9] = " << y0[9] << endl;		

		double end_value = y0[0] + 0.1;
		
		int counter = 0;
		while ( y0[0] < end_value ) 
		{
     		fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);
     			
     		if ( fehler != 0 ) 
			{
     			printf(" Gear4: error n° %d\n", 10 + fehler);
				break;
     		}
    
			cout << "t0: " << t0 << "; r: " << y0[0] << endl;

     		xend = sampling_time * (counter + 2);
     		aufrufe = 0;  // actual number of calls
  				
			counter++;
		}
	
		// freeing dipole memory
		//delete [] dipole;

		// length of dipole vector = number of samples
		//int npoints = ddipx.size(); 
		//int freqs_size = npoints / 2.0;		

		//if ( freqs_size != 0 )
		//{
			// given the number of points and sampling time we can calculate freqs vector
			//vector<double> freqs = linspace( 0.0, 1.0 / ( 2.0 * sampling_time ), freqs_size );
		
			// deleting first frequency == 0
			//freqs.erase( freqs.begin() );

			// due to 2pi inside Fourier transofrm
			//multiply_vector( freqs, 2 * M_PI ); 
			// transforming reverse atomic time units to cm^-1
			//multiply_vector( freqs, constants::HZTOCM / constants::ATU );

			//cout << ">> Processing " << ics[0] << " trajectory. npoints = " << npoints << endl;

			//vector<double> ddipx_out = fft( ddipx );
			//vector<double> ddipy_out = fft( ddipy );
			//vector<double> ddipz_out = fft( ddipz );	
	
			// auxiliary variables to store interim variables
			//double power;
			//vector<double> intensities; 

	  		//for ( int i = 0; i < freqs_size; i++ )
	  		//{
				//power = ddipx_out[i] + ddipy_out[i] + ddipz_out[i];
				//power = power / pow( freqs[i], 2 ) * exp_hkt ;	
				
				//intensities.push_back( power );	
			//}
		
			//show( "freqs", "ints", freqs, intensities );

			//MPI_Send( &freqs_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
			////cout << ">> Process " << world_rank << " sends package size = " << freqs_size << " to root." << endl;

			//MPI_Send( &freqs[0], freqs_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			////cout << ">> Process " << world_rank << " sends frequency package to root" << endl;

			//MPI_Send( &intensities[0], freqs_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
			////cout << ">> Process " << world_rank << " sends intensities package to root" << endl;
		//}
		//else
		//{
			//MPI_Send( &freqs_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
			////cout << ">> Process " << world_rank << " sends package size = " << freqs_size << " to root." << endl;
		//}
	}
}

int main( int argc, char* argv[] )
{
	// Initialize the MPI environment
	MPI_Init( &argc, &argv );

	// getting id of the current process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// getting number of running processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
	
	if ( world_rank == 0 ) 
	{
		master_code( world_size );
	}	
	else
	{
		slave_code( world_rank );
	}

	MPI_Finalize();
	
	return 0;
}
