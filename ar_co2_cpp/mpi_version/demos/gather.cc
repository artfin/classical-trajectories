#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>

// creates an array of random numbers; each has a value from 0 - 1
float *create_rand_nums(int num_elements)
{
	float *rand_nums = (float*) malloc(sizeof(float) * num_elements);
	assert(rand_nums != NULL);

	for ( int i = 0; i < num_elements; i++ )
	{
		rand_nums[i] = (rand() / (float) RAND_MAX);
	}
	return rand_nums;
}

int main(int argc, char** argv)
{
  // Seed the random number generator to get different results each time
  srand(time(NULL));

  MPI_Init(NULL, NULL);
  
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  const int LENGTH = 24;
  const int num_elements_per_proc = 6;

  // Create a random array of elements on the root process.
  float *rand_nums = NULL;
  if (world_rank == 0) {
      rand_nums = create_rand_nums(LENGTH);
  }

  // For each process, create a buffer that will hold a subset of the entire array
  float *sub_rand_nums = (float*) malloc(sizeof(float) * num_elements_per_proc);
  assert(sub_rand_nums != NULL);

  // Scatter the random numbers from the root process to all processes in the MPI WORLD
  MPI_Scatter(rand_nums, num_elements_per_proc, MPI_FLOAT, sub_rand_nums,
	      num_elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);
  // rand_nums -- array to be scattered
  // num_elements_per_proc -- number of elements to be distributed per process
  // MPI_FLOAT -- type of distributed element
  // sub_rand_nums -- name of the data container in the namespace of the process
  // 0 -- number of process to scatter data
  
  printf("Rank: %d; given numbers: %.4f %.4f %.4f %.4f %.4f %.4f\n", world_rank, sub_rand_nums[0], sub_rand_nums[1], sub_rand_nums[2], sub_rand_nums[3], sub_rand_nums[4], sub_rand_nums[5]);

  // Clean up
  if ( world_rank == 0 )
  {
      free(rand_nums);
  }

  free(sub_rand_nums);
  
  // clean closing of MPI 
  MPI_Finalize();
}
