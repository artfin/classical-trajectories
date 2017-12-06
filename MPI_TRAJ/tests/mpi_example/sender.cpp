#include "sender.h"

void Sender::send( void )
{
	MPI_Send( &x, 1, MPI_INT, 1, 0, MPI_COMM_WORLD );
	std::cout << "Sended a message to process 1" << std::endl;	
}


