#pragma once

#include <fftw3.h>
#include "aux.h"

#define REALPART 0
#define IMAGPART 1

std::vector<double> fft( std::vector<double> signal );
