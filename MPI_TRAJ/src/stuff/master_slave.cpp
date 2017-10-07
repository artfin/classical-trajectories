#include <stdio.h>
#include <iostream>
#include "mpi.h"

#include <fstream>

#include <string>
#include <sstream>

#define ASSIGNMENT_TAG  1 
#define RESULTS_TAG     2

#define RESULTS_SIZE    1
#define NUM_CATEGORIES  10

const string OUTPUT_DIR = "output/dips/"; 

using namespace std;

void analyze_category( int cat, double* results )
{
	*results = 1;
}

void save_results( double cat, double* results )
{
	cout << "Saving results: " << &results << endl;
}

int main( int argc, char* argv[] )
{
	MPI_Init( &argc, &argv );

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	
	int num_proc;
	MPI_Comm_size( MPI_COMM_WORLD, &num_proc );

	int namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name( processor_name, &namelen );

	cout << "Program starting on rank " << rank << " of " << num_proc << " on " << processor_name << endl;

	int num_cat = NUM_CATEGORIES;

	int p, proc;
	double cat, results, cats_done;

	MPI_Status status;

	int FILENUMBER_MIN = 1;
	int FILENUMBER_MAX = 5;

	// SLAVE PROCESS CODE
	if ( rank != 0 )
	{
		for ( cat = 0; cat != -1; )
		{
			MPI_Recv( &cat, 1, MPI_INT, 0, ASSIGNMENT_TAG, MPI_COMM_WORLD, &status);

			if ( cat != -1 )
			{
				cout << "rank: " << rank << "; received assinment of category " << cat << endl;
				analyze_category( cat, &results );
				cout << "rank: " << rank << "; sending results from category " << cat << " to master." << endl;

				MPI_Send( &results, RESULTS_SIZE, MPI_DOUBLE, 0, RESULTS_TAG, MPI_COMM_WORLD );
			}
			else 
			{
				cout << "rank: " << rank << "; received die signal from master" << endl;
			}
		}
	}
	// MASTER PROCESS CODE
	else
	{	
		int curr_filenumber = FILENUMBER_MAX;

		// initial assignments
		for ( int slave = 1; slave < num_proc; slave++, curr_filenumber++ )
		{
			MPI_Send( &curr_filenumber, 1, MPI_INT, p, ASSIGNMENT_TAG, MPI_COMM_WORLD );
			cout << "rank: " << rank << "; assigning file " << curr_filenumber << " to process " << slave << endl;
		}

		// obtain results & continue doling out work until exhausted
		for ( cats_done = 0; cats_done < num_cat; )
		{
			MPI_Recv( &results, RESULTS_SIZE, MPI_DOUBLE, MPI_ANY_SOURCE, RESULTS_TAG, MPI_COMM_WORLD, &status );
			proc = status.MPI_SOURCE;

			save_results( cat_assign[proc], &results );
			cats_done++;

			MPI_Send( &cat, 1, MPI_INT, proc, ASSIGNMENT_TAG, MPI_COMM_WORLD );

			cout << "rank: " << rank << "; assigning category " << cat << " to process " << proc << endl;

			cat_assign[proc] = cat;

			if ( cat != -1 )
			{
				cat++;
				if ( cat == num_cat ) cat = -1;
			}
		}
	}

	cout << "rank: " << rank << "; calling MPI_Finalize() " << endl;

	MPI_Finalize();

	return 0;
}
