#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

const int EXIT_TAG = 42; // exiting tag
const int ICPERTRAJ = 2; // number of initial conditions per trajectory
const int NTRAJ = 6; // number of trajectories to be calculated

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

	// on the root process
	if ( world_rank == 0 ) {
		printf("World size: %d\n", world_size);

		// initializing array of initial conditions and filling it with random numbers
		double *ics = new double[NTRAJ * ICPERTRAJ];
		fill_rand_nums(ics, NTRAJ * ICPERTRAJ);		
		
		// MAKE SURE THAT THAT ICS DIVIDES EQUALLY BETWEEN PROCESSES !
		
		int alive = world_size - 1;

		int source;
		int output;
		MPI_Status status;

		// sending initial conditions to slave processes
		int index = 0; // number of element in array to be sent 

		// sending first message to slaves
		for ( int i = 1; i < world_size; i++ ) {	
			
			MPI_Send(&ics[index], ICPERTRAJ, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			printf("Master process send message %.4f ... to %d\n", ics[index], 1);
			index += ICPERTRAJ; // sending data by chunks of ICPERTRAJ
		}

		while ( true ) {	
			
			if ( alive == 0 ) {
				break;
			}
			
			MPI_Recv(&output, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			source = status.MPI_SOURCE;
			printf("Master received report from %d\n", source);

			// sending new chunk of work
			if ( index < ICPERTRAJ * NTRAJ ) {

				MPI_Send(&ics[index], ICPERTRAJ, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);
				index += ICPERTRAJ;

			} else {

				// work is done, sending a killing message
				MPI_Send(&ics[index], ICPERTRAJ, MPI_DOUBLE, source, EXIT_TAG, MPI_COMM_WORLD);
				// decrementing alive processes counter
				alive--;
			}	
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
			
			// sleeping for a couple of seconds
			usleep(2000000);
			
			// sending report to master
			int report = 1;
			MPI_Send(&report, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
	
	MPI_Finalize();
}



