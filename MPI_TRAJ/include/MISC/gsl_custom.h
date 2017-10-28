#pragma once

#include <fstream>
#include <string>
#include <gsl/gsl_histogram.h>

void save_histogram( gsl_histogram *histogram, const int NBINS, const std::string filename  );
