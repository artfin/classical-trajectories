#pragma once

#include "aux.h"

#include <fftw3.h>
#include <algorithm>

#define REALPART 0
#define IMAGPART 1

std::vector<double> fft( std::vector<double> signal );
std::vector<double> fft_full( std::vector<double> &signal );

void fftfreq( std::vector<double> &freqs, const int len, const double d );
void fftshift( std::vector<double> &arr, const int len );
void ifftshift( std::vector<double> &arr, const int len );
