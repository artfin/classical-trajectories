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

const double LIGHTSPEED_CM = 3.0 * 1E10; // cm / s
const long double PLANCKCONST = 6.62 * 1E-34; // j * s
const long double PLANCKCONST_REDUCED = 1.054571800 *1E-34; // j * s
const long double LOSHMIDT_CONSTANT_CM = 2.6867774 * 1E19; // cm^-3
const double HZTOCM = 3.335631 * pow( 10, -11 ); // hz to cm
const double CMTOHZ = LIGHTSPEED_CM; // 
const long double BOLTZCONST = 1.23 * 1E-23; // j / k
const double Temperature = 295.0; // k

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

void plot_signal( Gnuplot &gp, vector< vector<double>> &freqs, vector< vector<double>> &v, vector< string > titles, bool set_logscale )
{
	//gp << "set terminal pngcairo enhanced font 'Dejavu Sans'\n";
	gp << "set encoding utf8\n";
	gp << "set format y '%.0sx10^{%S}'\n";
	gp << "set xrange [0.0 : 500];\n";
	gp << "set yrange [4.76e-79:4.76e-77];\n";	

	if ( set_logscale )
	{
		gp << "set logscale y;\n";
	}

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

vector<double> calculate_alpha( vector<double> &freqs, vector<double> &intensities )
{
	double alpha;
	vector<double> alphas;

	long double constant = 4 * pow(M_PI, 2) * pow(LOSHMIDT_CONSTANT_CM, 2) / (3 * PLANCKCONST_REDUCED * LIGHTSPEED_CM);

	long double in_braces;
	for ( int i = 0; i < freqs.size(); i++ )
	{
		in_braces = 1 - exp( - (PLANCKCONST_REDUCED * freqs[i] * CMTOHZ) / (BOLTZCONST * Temperature) );
		alpha = constant * freqs[i] * CMTOHZ * in_braces * intensities[i];	
		alphas.push_back( alpha );
	}

	return alphas;
}

int main( int argc, char* argv[] )
{
	bool make_specfunc = true;
	bool make_binary_absorption = false;

	if ( argc == 2 )
	{
		string str = argv[1];
		if ( str == "alpha" )
		{
			make_specfunc = false;
			make_binary_absorption = true;
		}
	}

	// ##############################################################
	vector<double> lbs, lbs2, lbs3, lbs4, lbs5, lbs6, lbs7, lbs8, lbs9, lbs10, lbs11, lbs12, lbs13, lbs14, lbs15, lbs16, lbs17;
	vector<double> ubs, ubs2, ubs3, ubs4, ubs5, ubs6, ubs7, ubs8, ubs9, ubs10, ubs11, ubs12, ubs13, ubs14, ubs15, ubs16, ubs17;
	vector<double> contents, contents2, contents3, contents4, contents5, contents6, contents7, contents8, contents9, contents10, contents11, contents12, contents13, contents14, contents15, contents16, contents17;

	vector<double> lbs_mcmc1, lbs_mcmc2, lbs_mcmc3, lbs_mcmc4, lbs_mcmc5;
	vector<double> ubs_mcmc1, ubs_mcmc2, ubs_mcmc3, ubs_mcmc4, ubs_mcmc5;
	vector<double> contents_mcmc1, contents_mcmc2, contents_mcmc3, contents_mcmc4, contents_mcmc5;

	vector<double> lbs_mcmc_gunsight1, lbs_mcmc_gunsight2, lbs_mcmc_gunsight3, lbs_mcmc_gunsight4, lbs_mcmc_gunsight5;
	vector<double> ubs_mcmc_gunsight1, ubs_mcmc_gunsight2, ubs_mcmc_gunsight3, ubs_mcmc_gunsight4, ubs_mcmc_gunsight5;
	vector<double> contents_mcmc_gunsight1, contents_mcmc_gunsight2, contents_mcmc_gunsight3, contents_mcmc_gunsight4, contents_mcmc_gunsight5;

	vector<double> lbs_mcmc_no_weight1, lbs_mcmc_no_weight2;
	vector<double> ubs_mcmc_no_weight1, ubs_mcmc_no_weight2;
	vector<double> contents_mcmc_no_weight1, contents_mcmc_no_weight2;

	vector<double> lbs_buryak1, lbs_buryak2, lbs_buryak3, lbs_buryak4, lbs_buryak5, lbs_buryak6, lbs_buryak7, lbs_buryak8, lbs_buryak9;
   	vector<double> ubs_buryak1, ubs_buryak2, ubs_buryak3, ubs_buryak4, ubs_buryak5, ubs_buryak6, ubs_buryak7, ubs_buryak8, ubs_buryak9;
	vector<double> contents_buryak1, contents_buryak2, contents_buryak3, contents_buryak4, contents_buryak5, contents_buryak6, contents_buryak7, contents_buryak8, contents_buryak9;	
	// ##############################################################

	// ##############################################################
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

	//read_file( "mcmc_gunsight_1_no_weight", lbs_mcmc_no_weight1, ubs_mcmc_no_weight1, contents_mcmc_no_weight1 );
	//read_file( "new_results/mcmc_gunsight_5_no_weight", lbs_mcmc_no_weight2, ubs_mcmc_no_weight2, contents_mcmc_no_weight2 );

	read_file( "new_results/buryak_4600_50_6_025_1000_2325_10", lbs_buryak1, ubs_buryak1, contents_buryak1 );
	read_file( "new_results/buryak_4600_25_6_025_1000_4625_10", lbs_buryak2, ubs_buryak2, contents_buryak2 );
	read_file( "new_results/buryak_4600_50_6_0125_1000_4557_10", lbs_buryak3, ubs_buryak3, contents_buryak3 );
	read_file( "new_results/buryak_4600_50_6_025_500_2325_10", lbs_buryak4, ubs_buryak4, contents_buryak4 );
	read_file( "new_results/buryak_4600_50_6_025_100_2325_10", lbs_buryak5, ubs_buryak5, contents_buryak5 );
	read_file( "new_results/buryak_4650_50_625_025_500_2444_10_simpson", lbs_buryak6, ubs_buryak6, contents_buryak6 );
	read_file( "new_results/buryak_4650_50_625_025_500_2444_10_trapezoid", lbs_buryak7, ubs_buryak7, contents_buryak7 );
	read_file( "new_results/buryak_4650_50_625_025_500_2444_20_simpson", lbs_buryak8, ubs_buryak8, contents_buryak8 );
	read_file( "new_results/buryak_9300_50_625_025_500_4862_20_simpson", lbs_buryak9, ubs_buryak9, contents_buryak9 );
	// ##############################################################
    
	// ##############################################################
	vector< vector<double>> x { lbs, lbs2, lbs3, lbs4, lbs5, lbs6, lbs7, lbs8, lbs9, lbs10, lbs11, lbs12, lbs13, lbs14, lbs15, lbs16, lbs17 };
	vector< vector<double>> y { contents, contents2, contents3, contents4, contents5, contents6, contents7, contents8, contents9, contents10, contents11, contents12, contents13, contents14, contents15, contents16, contents17 };

	vector< vector<double>> x_mcmc { lbs_mcmc2, lbs_mcmc3, lbs_mcmc4, lbs_mcmc5 };
	vector< vector<double>> y_mcmc { contents_mcmc2, contents_mcmc3, contents_mcmc4, contents_mcmc5 };

	//vector< vector<double>> x_mcmc_gunsight{ lbs_mcmc_gunsight2, lbs_mcmc_gunsight3, lbs_mcmc_gunsight4, lbs_mcmc_gunsight5, lbs_mcmc_no_weight1, lbs_mcmc_no_weight2 };
	//vector< vector<double>> y_mcmc_gunsight{ contents_mcmc_gunsight2, contents_mcmc_gunsight3, contents_mcmc_gunsight4, contents_mcmc_gunsight5, contents_mcmc_no_weight1, contents_mcmc_no_weight2 };

	vector<vector<double>> x_mcmc_gunsight{ lbs_mcmc_no_weight1 };
	vector<vector<double>> y_mcmc_gunsight{ contents_mcmc_no_weight1 };
	
	vector<vector<double>> x_buryak{ lbs_buryak1, lbs_buryak2, lbs_buryak3, lbs_buryak4, lbs_buryak5, lbs_buryak6, lbs_buryak7, lbs_buryak8, lbs_buryak9 };
	vector<vector<double>> y_buryak{ contents_buryak1, contents_buryak2, contents_buryak3, contents_buryak4, contents_buryak5, contents_buryak6, contents_buryak7, contents_buryak8, contents_buryak9 };

    // ##############################################################
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
    // ##############################################################

    // ##############################################################
	vector< string > titles_mcmc;
	//titles_mcmc.push_back( "MCMC ICS; 1.000 trajs (R0 = 40.0 a0)" );
	//titles_mcmc.push_back( "MCMC ICS; 50.000 trajs (R0 = 40.0 a0)" );
	//titles_mcmc.push_back( "MCMC ICS; 50.000 trajs (R0 = 200.0 a0)" );	
	//titles_mcmc.push_back( "MCMC ICS; 200.000 trajs (R0 = 200.0 a0)" );	
	//titles_mcmc.push_back( "MCMC ICS; 500.000 trajs (R0 = 200.0 a0)" );	
    // ##############################################################

    // ##############################################################
	vector< string > titles_mcmc_gunsight;
	//titles_mcmc_gunsight.push_back( "MCMC ICS; 1.000 trajs (b < 1.0, R0 = 200 a0)" );
	//titles_mcmc_gunsight.push_back( "MCMC ICS; 1.000 trajs (b < 2.0, R0 = 200 a0)" );
	//titles_mcmc_gunsight.push_back( "MCMC ICS; 1.000 trajs(b < 1.0, R0 = 50a0)" );
	//titles_mcmc_gunsight.push_back( "MCMC ICS; 5.000 trajs(b < 1.0, R0 = 50a0)" );
	//titles_mcmc_gunsight.push_back( "MCMC ICS; 3.000 trajs(b < 1.0, R0 = 200a0); NO WEIGHT" );
	//titles_mcmc_gunsight.push_back( "MCMC ICS; 1.000 trajs(b < 1.0, R0 = 50a0); NO WEIGHT" );
    // ##############################################################

    // ##############################################################
	vector<string> titles_buryak;
	titles_buryak.push_back( "VMAX: 4600; VSTEP: 50; BMAX: 6; BSTEP: 0.25; SAMPLING TIME: 1000, NTRAJS: 2325; BIN: 10" );
	titles_buryak.push_back( "VMAX: 4600; VSTEP: 25; BMAX: 6; BSTEP: 0.25; SAMPLING TIME: 1000, NTRAJS: 4625; BIN: 10" );
	titles_buryak.push_back( "VMAX: 4600; VSTEP: 50; BMAX: 6; BSTEP: 0.125; SAMPLING TIME: 1000, NTRAJS: 4557; BIN: 10" );
	titles_buryak.push_back( "VMAX: 4600; VSTEP: 50; BMAX: 6; BSTEP: 0.25; SAMPLING TIME: 500, NTRAJS: 2325; BIN: 10" );
	titles_buryak.push_back( "VMAX: 4600; VSTEP: 50; BMAX: 6; BSTEP: 0.25; SAMPLING TIME: 100, NTRAJS: 2325; BIN: 10" );
	titles_buryak.push_back( "VMAX: 4650; VSTEP: 50; BMAX: 6.25; BSTEP: 0.25; SAMPLING TIME: 500, NTRAJS: 2444; BIN: 10; SIMPSON" );
	titles_buryak.push_back( "VMAX: 4650; VSTEP: 50; BMAX: 6.25; BSTEP: 0.25; SAMPLING TIME: 500, NTRAJS: 2444; BIN: 10; TRAPEZOID" );
	titles_buryak.push_back( "VMAX: 4650; VSTEP: 50; BMAX: 6.25; BSTEP: 0.25; SAMPLING TIME: 500, NTRAJS: 2444; BIN: 20; SIMPSON" );
	titles_buryak.push_back( "VMAX: 9300; VSTEP: 50; BMAX: 6.25; BSTEP: 0.25; SAMPLING TIME: 500, NTRAJS: 4862; BIN: 20; SIMPSON" );
    // ##############################################################

	Gnuplot gp;
	//plot_signal( gp, x, y, titles );
	//plot_signal( gp, x_mcmc, y_mcmc, titles_mcmc );
	//plot_signal( gp, x_mcmc_gunsight, y_mcmc_gunsight, titles_mcmc_gunsight );

	if ( make_specfunc )
	{
		cout << "Plotting spectral function..." << endl;
		plot_signal( gp, x_buryak, y_buryak, titles_buryak, true );
	}
	
	if ( make_binary_absorption )
	{
		cout << "Plotting alpha..." << endl;

		vector<double> temp;
		vector<vector<double>> alpha_buryak;
	
		for ( int i = 0; i < x_buryak.size(); i++ )
		{
			temp = calculate_alpha( x_buryak[i], y_buryak[i] );
			alpha_buryak.push_back( temp );
			temp.clear();
		}

		plot_signal( gp, x_buryak, alpha_buryak, titles_buryak, false ); 
	}
	return 0;
}

