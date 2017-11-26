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
using std::get;
using std::make_pair;
using std::string;

void sample_sin( vector<double> &input, const double sampling_time, const int N )
{
	double x0 = 0.0;
	double x;	

	for ( int i = 0; i < N; i++ )
	{
		x = x0 + i * sampling_time; 
		
		input.push_back( 7 * sin( 2 * M_PI * x ) + 3 * sin( 2 * M_PI * 5 * x ) + 5 * sin(2 * M_PI * x * 0.5) + sin( 2 * M_PI * x * 0.1 ) ); 
	}
}

void show_vector( vector<double> &v, string name ) 
{
	for ( int i = 0; i < v.size(); i++ )
	{
		cout << name << "[" << i << "] = " << v[i] << endl;
	}
}	

void plot_signal( Gnuplot &gp, vector< vector<double>> &x, vector< vector<double>> &y, vector<string> &titles )
{
	gp << "set xrange [-2:2];\n";

	string gnuplot_cmd = "plot ";
	string temp;

	for ( int counter = 0; counter < titles.size(); counter++ )
	{
		temp = "'-' with lines title '" + titles[counter] + "'";

		if ( counter < titles.size() - 1 )
		{
			temp += ","; 
		}

		gnuplot_cmd += temp;
	}

	cout << "gnuplot_cmd: " << gnuplot_cmd << endl;

	gp << gnuplot_cmd << endl;

	vector< pair<double, double>> signal;

	for ( int i = 0; i < x.size(); i++ )
	{
		for ( int j = 0; j < x[i].size(); j++ )
		{
			signal.push_back( make_pair( x[i][j], y[i][j] ));
		}

		for ( int i = 0; i < signal.size(); i++ )
		{
			cout << "signal[i]: " << get<0>( signal[i] ) << " " <<
									 get<1>( signal[i] ) << endl;
		}
		
		gp.send1d( signal );

		signal.clear();
	}
	
	gp.flush();
}


int main( int argc, char* argv[] )
{
	vector<double> input;
	vector<double> freqs_one_side, output_one_side;
	vector<double> freqs_two_side, output_two_side;

	int N = 5000;
	double sampling_time = 0.1;	

	int freq_points_one_side = (int) (N + 1) / 2.0; 

	sample_sin( input, sampling_time, N );

	// ###########################################################
	fft_positive( input );
	// ###########################################################

	// ###########################################################
	freqs_one_side = linspace( 0.0, 1.0 / ( 2.0 * sampling_time ), freq_points_one_side );

	output_one_side = fft_one_side( input );

	multiply_vector( output_one_side, (double) 2.0 / N );
	// ###########################################################

	// ###########################################################
	fftfreq( freqs_two_side, N, sampling_time ); 
   	// freqs should IFFTSHIFTed! 
	ifftshift( freqs_two_side, N );	
	
	output_two_side = fft_two_side( input );
	// output should be IFFTSHIFT!
	ifftshift( output_two_side, N );
	
	multiply_vector( output_two_side, (double) 2.0 / N );
	// ###########################################################
	
	//show_vector( output_one_side, "output-one-side" );

	cout << "###################" << endl;

	//show_vector( freqs, "freqs" );

	vector< vector<double> > x{  freqs_one_side, freqs_two_side };
	vector< vector<double> > y{  output_one_side, output_two_side };
	vector< string > titles{ "one-side", "two-side" };

	Gnuplot gp;
	plot_signal( gp, x, y, titles ); 	

	return 0;
}
