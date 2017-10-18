#include <iostream>

#include <vector>
#include <algorithm>

#include <fstream>
#include <sstream>

using namespace std;

double median( vector<double> v )
{
	sort( v.begin(), v.end() );

	double med;	
	if ( v.size() % 2 == 0 )
	{
		med = 0.5 * ( v[ v.size() / 2.0 - 1 ] + 
			  		  v[ v.size() / 2.0 ]);
	}
	else
	{
		med = v[ v.size() / 2 ];	
	}

	return med;
}

vector<double> smooth( vector<double> contents )
{
	vector<double> slice;
	vector<double> filtered;

	double res;

	for ( int i = 0; i < contents.size(); i++ )
	{
		if ( i == 0 )
		{
			slice.push_back( contents[i] );
			slice.push_back( contents[i] );
			slice.push_back( contents[i + 1] );
		}
	
		else if ( i == contents.size() - 1 )
		{
			slice.push_back( contents[i - 1] );
			slice.push_back( contents[i] );
			slice.push_back( contents[i] );
		}
		else
		{
			slice.push_back( contents[i - 1] );
			slice.push_back( contents[i] );
			slice.push_back( contents[i + 1] );
		}

		res = median( slice );

		filtered.push_back( res );

		slice.clear();
	}

	return filtered;
}

void save_smoothed( vector<double> lb, vector<double> ub, vector<double> v )
{
	ofstream output( "spectrum_smoothed.txt" );
	
	if ( output.is_open() )
	{
		for ( int i = 0; i < v.size(); i++ )
		{
			output << lb[i] << " " << ub[i] << " " << v[i] << endl;
		}
	}

	output.close();
}

int main( int argc, char* argv[] )
{
	if ( argc != 2 )
	{
		cout << "USAGE: ./... (string) input filename" << endl;
		exit( 1 );
	}

	string input_filename = argv[1];

	ifstream input( input_filename.c_str() );

	string line;
	vector<double> v;
	double temp;

	vector<double> lower_bounds;
	vector<double> upper_bounds;
	vector<double> contents;

	if ( input.is_open() )
	{
		while ( getline( input, line ))
		{
			stringstream iss( line );
			
			while ( iss >> temp )
			{
				v.push_back( temp );
			}

			lower_bounds.push_back( v[0] );
			upper_bounds.push_back( v[1] );
			contents.push_back( v[2] );

			v.clear();
		}
		
		input.close();
	}
	else
	{
		cout << "Unable to open file." << endl;
	}

	contents = smooth( contents );	

	save_smoothed( lower_bounds, upper_bounds, contents );

	return 0;
}
