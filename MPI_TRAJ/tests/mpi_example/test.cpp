#include "mpi.h"
#include <iostream>

#include "sender.h"
#include "receiver.h"

using namespace std;

int main( int argc, char* argv[] )
{
	MPI_Init( &argc, &argv );

	int world_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );

	int world_size;
	MPI_Comm_size( MPI_COMM_WORLD, &world_size );

	if ( world_rank == 0 )
	{
		cout << "initializing Sender class" << endl;
		Sender sender( 5 );

		sender.send();
	}

	if ( world_rank == 1 )
	{
		cout << "initializing Receiver class" << endl;
		Receiver receiver;
		receiver.receive();
	}

	MPI_Finalize();

	return 0;
}
