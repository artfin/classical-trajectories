#include "receiver.h"

void Receiver::receive( void )
{
	MPI_Status status;
	MPI_Recv( &x, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status );
	std::cout << "Received message from process " << status.MPI_SOURCE << 
			": " << x << std::endl;
}


