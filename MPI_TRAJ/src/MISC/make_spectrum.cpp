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

void plot_signal( Gnuplot &gp, vector< vector<double>> &freqs, vector< vector<double>> &v, vector< string > titles )
{
	gp << "set xrange [0.0 : 500];\n";
	gp << "set yrange [0.005 : 5];\n";	

	gp << "set logscale y;\n";
		
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

	for ( int i = 0; i < freqs.size(); i++ )
	{
		for ( int j = 0; j < freqs[i].size(); j++ )
		{
			signal.push_back( make_pair( freqs[i][j], v[i][j] ));
		}

		gp.send1d( signal );

		signal.clear();
	}
	
	gp.flush();
}

int main()
{
	vector<double> lbs, lbs2, lbs3, lbs4, lbs5, lbs6, lbs7, lbs8, lbs9, lbs10, lbs11, lbs12, lbs13, lbs14, lbs15, lbs16, lbs17;
	vector<double> ubs, ubs2, ubs3, ubs4, ubs5, ubs6, ubs7, ubs8, ubs9, ubs10, ubs11, ubs12, ubs13, ubs14, ubs15, ubs16, ubs17;
	vector<double> contents, contents2, contents3, contents4, contents5, contents6, contents7, contents8, contents9, contents10, contents11, contents12, contents13, contents14, contents15, contents16, contents17;

	vector<double> lbs_mcmc1, lbs_mcmc2, lbs_mcmc3, lbs_mcmc4, lbs_mcmc5;
	vector<double> ubs_mcmc1, ubs_mcmc2, ubs_mcmc3, ubs_mcmc4, ubs_mcmc5;
	vector<double> contents_mcmc1, contents_mcmc2, contents_mcmc3, contents_mcmc4, contents_mcmc5;

	vector<double> lbs_mcmc_gunsight1, lbs_mcmc_gunsight2, lbs_mcmc_gunsight3, lbs_mcmc_gunsight4, lbs_mcmc_gunsight5;
	vector<double> ubs_mcmc_gunsight1, ubs_mcmc_gunsight2, ubs_mcmc_gunsight3, ubs_mcmc_gunsight4, ubs_mcmc_gunsight5;
	vector<double> contents_mcmc_gunsight1, contents_mcmc_gunsight2, contents_mcmc_gunsight3, contents_mcmc_gunsight4, contents_mcmc_gunsight5;

	read_file( "new_results/weights_test", lbs, ubs, contents );
	read_file( "new_results/weights_test2", lbs2, ubs2, contents2 );
	read_file( "new_results/weights_test3", lbs3, ubs3, contents3 );
	read_file( "new_results/weights_test4", lbs4, ubs4, contents4 );
	read_file( "new_results/weights_test5", lbs5, ubs5, contents5 );
	read_file( "new_results/weights_test6", lbs6, ubs6, contents6 );
	read_file( "new_results/weights_test7", lbs7, ubs7, contents7 );
	read_file( "new_results/weights_test8", lbs8, ubs8, contents8 );
	read_file( "new_results/weights_test9", lbs9, ubs9, contents9 );
	read_file( "new_results/weights_test10", lbs10, ubs10, contents10 );
	read_file( "new_results/weights_test11", lbs11, ubs11, contents11 );
	read_file( "new_results/weights_test12", lbs12, ubs12, contents12 );
	read_file( "new_results/weights_test13", lbs13, ubs13, contents13 );
	read_file( "new_results/weights_test14", lbs14, ubs14, contents14 );
	read_file( "new_results/weights_test15", lbs15, ubs15, contents15 );
	read_file( "new_results/weights_test16", lbs16, ubs16, contents16 );
	read_file( "new_results/weights_test17", lbs17, ubs17, contents17 );

	read_file( "new_results/mcmc1", lbs_mcmc1, ubs_mcmc1, contents_mcmc1 );
	read_file( "new_results/mcmc2", lbs_mcmc2, ubs_mcmc2, contents_mcmc2 );
	read_file( "new_results/mcmc3", lbs_mcmc3, ubs_mcmc3, contents_mcmc3 );
	read_file( "new_results/mcmc4", lbs_mcmc4, ubs_mcmc4, contents_mcmc4 );
	read_file( "new_results/mcmc5", lbs_mcmc5, ubs_mcmc5, contents_mcmc5 );

	read_file( "new_results/mcmc_gunsight_1", lbs_mcmc_gunsight1, ubs_mcmc_gunsight1, contents_mcmc_gunsight1 );
	read_file( "new_results/mcmc_gunsight_2", lbs_mcmc_gunsight2, ubs_mcmc_gunsight2, contents_mcmc_gunsight2 );
	read_file( "new_results/mcmc_gunsight_3", lbs_mcmc_gunsight3, ubs_mcmc_gunsight3, contents_mcmc_gunsight3 );
	read_file( "new_results/mcmc_gunsight_4", lbs_mcmc_gunsight4, ubs_mcmc_gunsight4, contents_mcmc_gunsight4 );
	read_file( "new_results/mcmc_gunsight_5", lbs_mcmc_gunsight5, ubs_mcmc_gunsight5, contents_mcmc_gunsight5 );

	vector< vector<double>> x { lbs, lbs2, lbs3, lbs4, lbs5, lbs6, lbs7, lbs8, lbs9, lbs10, lbs11, lbs12, lbs13, lbs14, lbs15, lbs16, lbs17 };
	vector< vector<double>> y { contents, contents2, contents3, contents4, contents5, contents6, contents7, contents8, contents9, contents10, contents11, contents12, contents13, contents14, contents15, contents16, contents17 };

	vector< vector<double>> x_mcmc { lbs_mcmc2, lbs_mcmc3, lbs_mcmc4, lbs_mcmc5 };
	vector< vector<double>> y_mcmc { contents_mcmc2, contents_mcmc3, contents_mcmc4, contents_mcmc5 };

	vector< vector<double>> x_mcmc_gunsight{ lbs_mcmc_gunsight2, lbs_mcmc_gunsight3, lbs_mcmc_gunsight4, lbs_mcmc_gunsight5 };
	vector< vector<double>> y_mcmc_gunsight{ contents_mcmc_gunsight2, contents_mcmc_gunsight3, contents_mcmc_gunsight4, contents_mcmc_gunsight5 };

	vector< string > titles;
	titles.push_back( "PR in [0, -0.5](step=0.05); PT in [0, 1.0](step=0.05); R0 = 40.0" );
	titles.push_back( "PR in [0, -1.0](step=0.05); PT in [0, 1.0](step=0.05); R0 = 40.0" );
	titles.push_back( "PR in [0, -2.0](step=0.05); PT in [0, 1.0](step=0.05); R0 = 40.0" );
	titles.push_back( "PR in [0, -5.0](step=0.05); PT in [0, 1.0](step=0.05); R0 = 40.0" );
	titles.push_back( "PR in [0, -5.0](step=0.01); PT in [0, 1.0](step=0.05); R0 = 40.0" );
	titles.push_back( "PR in [0, -5.0](step=0.05); PT in [0, 5.0](step=0.05); R0 = 40.0" );
	titles.push_back( "PR in [0, -5.0](step=0.05); PT in [0, 1.0](step=0.05); R0 = 60.0" );
	titles.push_back( "PR in [0, -5.0](step=0.05); PT in [0, 1.0](step=0.05); R0 = 80.0" );
	titles.push_back( "PR in [0, -5.0](step=0.05); PT in [0, 1.0](step=0.05); R0 = 100.0" );
	//titles.push_back( "PR in [0, -10.0](step=0.05); PT in [0, 1.0](step=0.05); R0 = 200.0" );
	//titles.push_back( "PR in [0, -15.0](step=0.05); PT in [0, 1.0](step=0.05); R0 = 200.0" );
	//titles.push_back( "PR in [-5.0, -5.0](step=0.05); PT in [0, 10.0](step=0.05); R0 = 200.0" );
	//titles.push_back( "PR in [-10.0, -10.0](step=0.05); PT in [0, 10.0](step=0.05); R0 = 200.0" );
	//titles.push_back( "PR in [-1.0, -1.0](step=0.05); PT in [0, 30.0](step=0.05); R0 = 200.0" );
	//titles.push_back( "PR in [0.0, -15.0](step=0.05); PT in [1.0, 2.0](step=0.05); R0 = 200.0" );
	//titles.push_back( "PR in [0.0, -15.0](step=0.05); PT in [2.0, 5.0](step=0.05); R0 = 200.0" );
	//titles.push_back( "PR in [0.0, -15.0](step=0.05); PT in [5.0, 15.0](step=0.05); R0 = 200.0" );

	vector< string > titles_mcmc;
	//titles_mcmc.push_back( "MCMC ICS; 1.000 trajs (R0 = 40.0 a0)" );
	titles_mcmc.push_back( "MCMC ICS; 50.000 trajs (R0 = 40.0 a0)" );
	titles_mcmc.push_back( "MCMC ICS; 50.000 trajs (R0 = 200.0 a0)" );	
	titles_mcmc.push_back( "MCMC ICS; 200.000 trajs (R0 = 200.0 a0)" );	
	titles_mcmc.push_back( "MCMC ICS; 500.000 trajs (R0 = 200.0 a0)" );	

	vector< string > titles_mcmc_gunsight;
	titles_mcmc_gunsight.push_back( "MCMC ICS; 1.000 trajs (b < 1.0, R0 = 200 a0)" );
	titles_mcmc_gunsight.push_back( "MCMC ICS; 1.000 trajs (b < 2.0, R0 = 200 a0)" );
	titles_mcmc_gunsight.push_back( "MCMC ICS; 1.000 trajs(b < 1.0, R0 = 50a0)" );
	titles_mcmc_gunsight.push_back( "MCMC ICS; 5.000 trajs(b < 1.0, R0 = 50a0)" );

	Gnuplot gp;
	//plot_signal( gp, x, y, titles );
	//plot_signal( gp, x_mcmc, y_mcmc, titles_mcmc );
	plot_signal( gp, x_mcmc_gunsight, y_mcmc_gunsight, titles_mcmc_gunsight );

	return 0;
}

