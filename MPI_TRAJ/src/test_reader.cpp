#include <iostream>
#include <fstream>

// to perform string conversion
#include <string>
#include <sstream>

using namespace std;

int main( int argc, char* argv[] )
{
	if ( argc != 2 )
	{
		cout << "USAGE: ./.. (int) number_of_dipole file to process" << endl;
		exit( 1 );
	}

	int n = atoi( argv[1] );

	ostringstream strs;
	strs << n;
	string filename = "output/dips/" + strs.str() + ".bin";

	ifstream file( filename, ios::in | ios::binary | ios::ate );

	streampos size = file.tellg();
	int ndoubles = size / 8;

	cout << "Size = " << size << endl;
	cout << "Number of doubes: " << ndoubles << endl;

	char* memblock = new char [size];
	file.seekg( 0, ios::beg );
	file.read( memblock, size );
	file.close();

	cout << "The entire file content is in memory" << endl;
	double* double_values = (double*) memblock;

	for ( int i = 0; i < ndoubles / 2 + 1; i++ )
	{
		cout << "freq: " << double_values[i] << "; i: " << double_values[i + 1] << endl;
		i++;
	}

	return 0;
}
