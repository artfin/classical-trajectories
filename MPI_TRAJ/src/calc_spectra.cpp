#include <iostream>
#include <fstream>

// string conversion
#include <string>
#include <sstream>

#include <vector>
#include <algorithm>

#include <chrono>

using namespace std;

bool compareFunc( pair<double, double> &a, pair<double, double> &b )
{
	return a.first > b.first;
}

int main( int argc, char* argv[] )
{
	if ( argc != 3 )
	{
		cout << "USAGE: ./.. (int) starting file number (int) ending file number" << endl;
		exit( 1 );
	}

	int start = atoi( argv[1] );
	int end = atoi( argv[2] );

	auto startTime = chrono::high_resolution_clock::now();

	vector<double> freqs;
	vector<double> intensities;

	for ( int n = start; n < end + 1; n++ )
	{
		ostringstream strs;
		
		strs << n;
		string filename = "output/dips/" + strs.str() + ".bin";

		//cout << "filename: " << filename << endl;

		ifstream file( filename, ios::in | ios::binary | ios::ate );

		streampos size = file.tellg();
		int ndoubles = size / 8;

		//cout << "Size = " << size << endl;
		//cout << "Number of doubes: " << ndoubles << endl;

		char* memblock = new char [size];
		file.seekg( 0, ios::beg );
		file.read( memblock, size );
		file.close();

		//cout << "The entire file content is in memory" << endl;
		double* arr = (double*) memblock;
		
		if ( ndoubles == 0 )
		{
				//cout << "File is empty. Continuing..." << endl;
			continue;
		}

		int counter = 0;

		// current frequency and current intensity
		double curr_freq, curr_int;

		while ( counter < ndoubles )
		{
			curr_freq = arr[counter];
			curr_int = arr[counter + 1];

			vector<double>::iterator iter = find( freqs.begin(), freqs.end(), curr_freq );
			size_t index = distance( freqs.begin(), iter );
		
			// means that element is found 	
			if ( index != freqs.size() )
			{
				intensities[ index ] += curr_int; 
			}
			// adding freq to freqs vector
			else
			{
				freqs.push_back( curr_freq );
				intensities.push_back( curr_int );
			}
				
			counter += 2;
		}
	}

	// creating a pair
	vector< pair<double, double> > data ( freqs.size() );
	for ( size_t i = 0; i < freqs.size(); i++ )
	{
		data[i] = make_pair( freqs[i], intensities[i] );
	}
	
	// sorting by first component
	sort( data.begin(), data.end(), compareFunc );

	for ( int i = 0; i < freqs.size(); i++ )
	{
		cout << data[i].first << " " << data[i].second << endl; 
	}

	auto endTime = chrono::high_resolution_clock::now();

	cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds> ( endTime - startTime ).count() / 1000.0 << "s" << endl;;

	return 0;
}
