#include <iostream>
#include <fstream>

// string conversion
#include <string>
#include <sstream>

#include <vector>

#include <algorithm>

using namespace std;

template <typename T, typename Compare>
vector<size_t> sort_permutation(
	const vector<T>& vec, Compare& compare )
{
	vector<size_t> p( vec.size() );
	iota( p.begin(), p.end(), 0 );
	sort ( 	p.begin(), p.end(),
		[&] ( size_t i, size_t j )
		{
			return compare( vec[i], vec[j] );
		}	
	);

	return p;
}

template <typename T>
vector<T> apply_permutation(
	const vector<T>& vec,
	const vector<size_t>& p)
{
	vector<T> sorted_vec( vec.size() );
	transform( 
		p.begin(), p.end(), sorted_vec.begin(),
		[&] ( size_t i )
		{
			return vec[i];
		}
	);

	return sorted_vec;
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

	vector<double> freqs;
	vector<double> intensities;

	for ( int n = start; n < end + 1; n++ )
	{
		ostringstream strs;
		
		strs << n;
		string filename = "output/dips/" + strs.str() + ".txt";

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

	auto p = sort_permutation( freqs,  );
	
	freqs = apply_permutation( freqs, p );
	intensities = apply_permutation( intensities, p );	

	for ( int i = 0; i < freqs.size(); i++ )
	{
			//cout << "i: " << i << "; freq: " << freqs[i] << "; intensity: " << intensities[i] << endl;
		cout << freqs[i] << " " << intensities[i] << endl;
	}

	return 0;
}
