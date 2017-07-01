#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const int ICPERTRAJ = 2; // number of initial conditions per trajectory
const int NTRAJ = 4; // number of trajectories to be calculated

void fill_rand_nums(double* arr, int length ) {
    
	srand(time(NULL));	
	for ( int i = 0; i < length; i++ ) {
		arr[i] = (rand() / (double) RAND_MAX);
		printf("arr[%d] = %.4f\n", i, arr[i]);
	}
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
	printf("World rank: %d\n", world_rank);

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
		
		// sending initial conditions to slave processes
		// i -- number of process we're sending ics to
		int index = 0; // number of element in array to be sent  
		for ( int i = 1; i < world_size; i++ ) {
			
			MPI_Send(&ics[index], ICPERTRAJ, MPI_DOUBLE, i, 0, MPI_COMM_WORLD); 
			index += ICPERTRAJ; // sending data in chunks of ICPERTRAJ
			
			printf("Root process send message %.4f ... to %d\n", ics[index], i);
		}
	} 

	// on the slave process
	else {
		// array of initial conditions for slave process
		double *ics_sub = new double [ICPERTRAJ];
		
		MPI_Recv(&ics_sub[0], ICPERTRAJ, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
		printf("Process %d received message from root %.4f... \n", world_rank, ics_sub[0]);

		print_array(ics_sub, ICPERTRAJ);
		
		// sending back that process has finished
		int status = 1;
		printf("Process %d sending root a message with status %d\n", world_rank, status); 
		MPI_Send(&status, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
}



