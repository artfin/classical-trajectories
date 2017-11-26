#pragma once

#include <fftw3.h>
#include <vector>

std::vector<double> linspace( const double min, const double max, const int npoints );
void multiply_vector( std::vector<double> &v, const double factor );
void copy_to( std::vector<double> &v, fftw_complex* arr );
void copy_to( std::vector<double> &v, double* arr );
