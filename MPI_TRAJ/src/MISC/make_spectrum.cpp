#include <iostream>
#include <vector>

#include <string>
#include <cmath>
#include <fstream>
#include <sstream>

// Warn about use of deprecated functions
#define GNUPLOT_DEPRECATE_WARN
#include "gnuplot-iostream.h"

using namespace std;

const double LIGHTSPEED = 3.0 * 1E10;
const long double PLANCKCONST = 6.62 * 1E-34;
const long double BOLTZCONST = 1.23 * 1E-23;
const double Temperature = 300.0;

void read_file( string filename, vector<double> &lbs, 
				vector<double> &ubs, vector<double> &contents )
{
	ifstream input( filename );
	
	string line;
	double temp;
	vector<double> v;

	if ( input.is_open() )
	{
		while ( getline( input, line ) )
		{
			stringstream iss( line );

			while ( iss >> temp )
			{
				v.push_back( temp );
			}

			lbs.push_back( v[0] );
			ubs.push_back( v[1] );
			contents.push_back( v[2] );

			v.clear();
		}
	
		input.close();
	} 
	else 
	{
		cout << "Unable to open the file." << endl;
	}
}

void show_vector( string name, vector<double> v )
{
	cout << "##################################" << endl;
	
	for ( int i = 0; i < v.size(); i++ )
	{
		cout << name << "[" << i << "] = " << v[i] << endl;
	}
	
	cout << "##################################" << endl;
}

void plot_signal( Gnuplot &gp, vector<double> &freqs, vector<double> &v )
{
	vector< pair<double, double>> signal;

	for ( int i = 0; i < freqs.size(); i++ )
	{
		signal.push_back( make_pair( freqs[i], v[i] ));
	}
	
	gp << "set xrange [0:400]\n;";
		
	gp << "plot '-' with lines title 'signal'" << endl;
	gp.send1d( signal );

	gp.flush();
}

int main()
{
	vector<double> lbs;
	vector<double> ubs;
	vector<double> contents;

	read_file( "test", lbs, ubs, contents );

	Gnuplot gp;
	plot_signal( gp, lbs, contents );

	return 0;
}

