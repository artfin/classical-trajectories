#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "../jean-pierre/basis.h"
#include "../jean-pierre/vmblock.h"
#include "../jean-pierre/gear.h"
#include "../jean-pierre/t_dgls.h"

#include "../arco2_spectrum/matrix.h"
#include <fftw3.h>

#include <vector>

// to perform double to string conversion
#include <sstream>
#include <string>

// macros for real and imaginary parts
#define REALPART 0
#define IMAGPART 1

const int EXIT_TAG = 42; // exiting tag
const int ICPERTRAJ = 8; // number of initial conditions per trajectory

using namespace std;

void syst(REAL t, REAL *y, REAL *f);

int main()
{
	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

	// getting id of the current process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// getting number of running processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
	
	// number of alive processes
	int alive = world_size - 1;

	// auxiliary variables
	int source;
	int output;
	MPI_Status status;
	
	// on the root process
	if ( world_rank == 0 ) 
	{	
		FILE* inputfile = fopen("ics.txt", "r");

		double R, Theta;
		double *ics = new double[ICPERTRAJ];

		// result of reading values from file
		int scanfResult;

		// sending first message to slaves
		for ( int i = 1; i < world_size; i++ ) 
		{
			scanfResult = fwscanf(inputfile, L"%lf %lf %lf %lf %lf %lf %lf %lf\n", &ics[0], &ics[1], &ics[2], &ics[3], &ics[4], &ics[5], &ics[6], &ics[7]);
			printf("MASTER reads IC number %.1lf : %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", ics[0], ics[1], ics[2], ics[3], ics[4], ics[5], ics[6], ics[7]);  

			MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
		
		while ( true )
		{	
			// exits when all slaves are killed
			if ( alive == 0 )
			{
				break;
			}

			// receiving message from any of slaves
			MPI_Recv(&output, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			source = status.MPI_SOURCE;
			printf("MASTER received report from %d\n", source);
			
			// reading another line from file
			scanfResult = fwscanf(inputfile, L"%lf %lf %lf %lf %lf %lf %lf %lf\n", &ics[0], &ics[1], &ics[2], &ics[3], &ics[4], &ics[5], &ics[6], &ics[7]);	
			// if it's not ended yet then sending new chunk of work to slave
			if ( scanfResult != -1)
	   		{
				printf("MASTER reads IC number %.1lf: %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", ics[0], ics[1], ics[2], ics[3], ics[4], ics[5], ics[6], ics[7]);
				MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);
			}
			// work is done, sending a killing message
			else
			{
				MPI_Send(&ics[0], ICPERTRAJ, MPI_DOUBLE, source, EXIT_TAG, MPI_COMM_WORLD);
				alive--;
			}
		}
		
		fclose(inputfile);
	}

	// on the slave process
	else
	{
		REAL epsabs;    // absolute error bound
		REAL epsrel;    // relative error bound    
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

			if ( status.MPI_TAG == EXIT_TAG )
			{
				printf("Process %d exits work loop.\n", world_rank);
				break;
			}
			
			// output files
			ostringstream strs;
			strs << ics[0];
			string traj_filepath = "trajs/" + strs.str() + ".txt";
			cout << "Traj filepath: " << traj_filepath;

			string dip_filepath = "dips/" + strs.str() + ".txt";
			cout << "Dip filepath: " << dip_filepath;
			
			// .c_str() converts string to const char*
			FILE *trajectory_file = fopen(traj_filepath.c_str(), "w");
			// FILE *dipfft = fopen(dip_filepath, "w");

			N = 6;
			vmblock = vminit();
			y0 = (REAL*) vmalloc(vmblock, VEKTOR, N, 0);
			
			// out of memory?
  			if (! vmcomplete(vmblock))
			{ 
    			printf("mgear: out of memory.\n");
    			return 0;
 			}

		    double step = 3000;
  			double end = 700000;

  			int num_points = (end / step);

  			epsabs = 1E-13;
  			epsrel = 1E-13;
  
  			t0 = 0.0;
			
  			h = 0.1;         // initial step size
  			xend = step;     // initial right bound of integration
  			fmax = 1000000;  // maximal number of calls 

  			double *dipole = new double [3];

  			vector<double> ddipx(num_points);
  			vector<double> ddipy(num_points);
  			vector<double> ddipz(num_points);
			
			y0[0] = ics[1];
			y0[1] = ics[2];
			y0[2] = ics[3];
			y0[3] = ics[4];
			y0[4] = ics[5];
			y0[5] = ics[6];
  			y0[6] = ics[7];

			for (int i = 0; i < num_points; ++i)
  			{
     			fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);
     
     			if ( fehler != 0 ) 
				{
     				printf(" Gear4: error n° %d\n", 10 + fehler);
					break;
     			}
     
     			hamiltonian(dipole, y0[0], y0[1], y0[2], y0[3], y0[4], y0[5], y0[6], true);
     
     			fprintf(trajectory_file, "%f %.12f %.12f %.12f %.12f %.12f\n", t0, y0[0], y0[1], dipole[0], dipole[1], dipole[2]);

     			ddipx.push_back(dipole[0]);
     			ddipy.push_back(dipole[1]);
     			ddipz.push_back(dipole[2]);

     			xend = step * (i + 2);
     			aufrufe = 0;  // actual number of calls
  			}
 
			fclose(trajectory_file);

			// sending a report to master
			int report = 1;
			MPI_Send(&report, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}

	MPI_Finalize();
	
	return 0;
}

void syst (REAL t, REAL *y, REAL *f)
{
  (void)(t); // avoid unused parameter warning 

  double *out = new double[6];
  rhs(out, y[0], y[1], y[2], y[3], y[4], y[5], y[6]);
  // R  Theta pR pT phi theta J

  f[0] = out[0]; // dR/dt  
  f[1] = out[1]; // d(Theta)/dt
  f[2] = out[2]; // d(pR)/dt
  f[3] = out[3]; // d(pT)/dt
  f[4] = out[4]; // d(phi)/dt
  f[5] = out[5]; // d(theta)/dt

  delete [] out;
}


