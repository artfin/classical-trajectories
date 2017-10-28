#include "gsl_custom.h"

void save_histogram( gsl_histogram *histogram, const int NBINS, const std::string filename )
{
	std::ofstream file( filename );

	double lower_bound, higher_bound, bin_content;
	for ( int counter = 0; counter < NBINS; counter++ )
	{
		gsl_histogram_get_range( histogram, counter, &lower_bound, &higher_bound );
		bin_content = gsl_histogram_get( histogram, counter );

		file << lower_bound << " " << higher_bound << " " << bin_content << std::endl;
	}
	
	file.close();	
}

