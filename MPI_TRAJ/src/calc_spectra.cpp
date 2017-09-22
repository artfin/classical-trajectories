#include <iostream>
#include <fstream>

#include <cmath>

// string conversion
#include <string>
#include <sstream>

#include <vector>
#include <algorithm>

#include <chrono>

#include "psp_pes.h"
#include <hamiltonian.hpp>

using namespace std;

// hartree to joules
const double HTOJ = 4.35974417 * pow(10, -18);
// boltzmann constant
const long double BOLTZCONST = 1.38064852 * pow(10, -23);

const double Temperature = 70; // K

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

	string line;
	ifstream myfile( "input/ics.txt" );

	vector<double> theta( 50000 );
	vector<double> pr( 50000 );
	vector<double> ptheta( 50000 );
	vector<double> jphi( 50000 );
	vector<double> jtheta( 50000 );
	vector<double> jtot( 50000 );

	vector<double> ham_values( 50000 );

	double RDIST = 20.0;
	double jx, jy, jz;

	int i = 0;

	if ( myfile.is_open() )
	{
		while ( getline(myfile, line) )
		{
			istringstream sin(line);

			double temp;
			
			// it is number of trajectory N and R
			sin >> temp;
			sin >> temp;
			
			sin >> temp;
			theta[i] = temp;

			sin >> temp;
		   	pr[i] = temp;

			sin >> temp; 
			ptheta[i] = temp;

			sin >> temp;
		   	jphi[i] = temp;

			sin >> temp;
			jtheta[i] = temp;

			sin >> temp;
			jtot[i] = temp;

			i++;
		}

		myfile.close();
	}
	else
	{
		cout << "Unable to open file" << endl;
	}

	double h;
	for ( int i = 0; i < theta.size(); i++ )
	{
		jx = jtot[i] * sin(jtheta[i]) * cos(jphi[i]);
		jy = jtot[i] * sin(jtheta[i]) * sin(jphi[i]);
		jz = jtot[i] * cos(jtheta[i]);

		h = ham_value(RDIST, theta[i], pr[i], ptheta[i], jx, jy, jz);
		ham_values[i] = h;
	}	

	auto startTime = chrono::high_resolution_clock::now();

	double exp_hkt;

	vector<double> freqs;
	vector<double> intensities;

	for ( int n = start; n < end + 1; n++ )
	{
		exp_hkt = exp( - ham_values[n] * HTOJ / (BOLTZCONST * Temperature) );
		
		ostringstream strs;
		
		strs << n;
		string filename = "first_exp/dips/" + strs.str() + ".bin";

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

			//vector<double>::iterator iter = find( freqs.begin(), freqs.end(), curr_freq );
			//size_t index = distance( freqs.begin(), iter );
		
			// means that element is found 	
			//if ( index != freqs.size() )
			//{
			//intensities[ index ] += curr_int; 
			//}
			// adding freq to freqs vector
			//else
			//{
			freqs.push_back( curr_freq );
			intensities.push_back( curr_int * exp_hkt );
			//}
				
			counter += 2;
		}
	}

	cerr << "Binary files are read." << endl;

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
