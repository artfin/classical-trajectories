#include <iostream>
#include <cmath>
#include <vector>
#include <string>

#include "fft.h"

// Warn about use of deprecated functions
#define GNUPLOT_DEPRECATE_WARN
#include "gnuplot-iostream.h"

using std::cout;
using std::endl;

using std::vector;
using std::pair;
using std::make_pair;
using std::string;

void sample_sin( vector<double> &time, vector<double> &input, const double sampling_time, const int N )
{
	double x0 = 0.0;
	double x;	

	for ( int i = 0; i < N; i++ )
	{
		x = x0 + i * sampling_time; 
		
		time.push_back( x );
		input.push_back( 7 * sin( 2 * M_PI * x ) + 3 * sin( 2 * M_PI * 5 * x ) ); 
	}
}

void plot_signal( Gnuplot &gp, vector<double> &freqs, vector<double> &output )
{
	gp << "plot '-' with lines title 'signal'\n"; 
	
	vector< pair<double, double>> signal;

	for ( int i = 0; i < freqs.size(); i++ )
	{
		signal.push_back( make_pair( freqs[i], output[i] ));
	}
	
	gp.send1d( signal );
	gp.flush();
}

void show_vector( vector<double> &v, string name ) 
{
	for ( int i = 0; i < v.size(); i++ )
	{
		cout << name << "[" << i << "] = " << v[i] << endl;
	}
}	

int main( int argc, char* argv[] )
{
	vector<double> time, input;
	vector<double> freqs, output;

	int N = 100;
	double sampling_time = 0.05;	
	
	sample_sin( time, input, sampling_time, N );
	
	//freqs = linspace( -1.0 / (2.0 * sampling_time), 1.0 / ( 2.0 * sampling_time ), N );

	output = fft_full( input );
	ifftshift( output, N );
	//multiply_vector( output, (double) 1.0 / pow(N, 2) / 2.0 );
	//show_vector( output, "output" );

	cout << "###################" << endl;

	fftfreq( freqs, N, sampling_time ); 
   	ifftshift( freqs, N );	
	//show_vector( freqs, "freqs" );

	Gnuplot gp;
	plot_signal( gp, freqs, output ); 	

	return 0;
}
