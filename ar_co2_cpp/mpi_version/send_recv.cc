#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const int EXIT_TAG = 42; // exiting tag
const int ICPERTRAJ = 2; // number of initial conditions per trajectory
const int NTRAJ = 3; // number of trajectories to be calculated

void fill_rand_nums(double* arr, int length ) {
    
	srand(time(NULL));	
	printf("-------------\n");
	for ( int i = 0; i < length; i++ ) {
		arr[i] = (rand() / (double) RAND_MAX);
		printf("arr[%d] = %.4f\n", i, arr[i]);
	}
	printf("-------------\n");
}

void print_array(double* arr, int length) {
	
	for ( int i = 0; i < length; i++ ) {
		printf("arr[%d] = %.4lf\n", i, arr[i]);
	}
}

int main()
{
	// Initialize the MPI environment
	MPI_Init(NULL, NULL);

	// getting id of the current process
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	// printf("World rank: %d\n", world_rank);

	// getting number of running processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size); 

	const int NCYCLES = (int) NTRAJ / (world_size - 1);

	// on the root process
	if ( world_rank == 0 ) {
		printf("World size: %d\n", world_size);
		printf("Number of cycles per slave: %d\n", NCYCLES);

		// initializing array of initial conditions and filling it with random numbers
		double *ics = new double[NTRAJ * ICPERTRAJ];
		fill_rand_nums(ics, NTRAJ * ICPERTRAJ);		
		
		// MAKE SURE THAT THAT ICS DIVIDES EQUALLY BETWEEN PROCESSES !

		// sending initial conditions to slave processes
		int index = 0; // number of element in array to be sent  
		while ( index < NTRAJ * ICPERTRAJ ) {
			
			// i -- number of process we're sending ics to
			for ( int i = 1; i < world_size; i++ ) {
				
					MPI_Send(&ics[index], ICPERTRAJ, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
					printf("Master process send message %.4f ... to %d\n", ics[index], i);
		    		index += ICPERTRAJ; // sending data by chunks of ICPERTRAJ
			}
		}
		
		// sending killing message
		for ( int i = 1; i < world_size; i++ ) {
			MPI_Send(&ics[index], ICPERTRAJ, MPI_DOUBLE, i, EXIT_TAG, MPI_COMM_WORLD);
		}
	} 

	// on the slave process
	else {
		// array of initial conditions for slave process
		double *ics_sub = new double [ICPERTRAJ];
		
		while (true) {

			MPI_Status status;
			MPI_Recv(&ics_sub[0], ICPERTRAJ, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			if ( status.MPI_TAG == EXIT_TAG ) {
				printf("Process %d exits work loop.\n", world_rank);
				break;
			}
	
			printf("Process %d received message from root %.4f... \n", world_rank, ics_sub[0]);
		}
	}

	MPI_Finalize();
}



