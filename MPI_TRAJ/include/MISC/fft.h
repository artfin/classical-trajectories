#pragma once

#include "aux.h"

#include <fftw3.h>
#include <algorithm>
#include <fstream>

#include <iostream>

#define REALPART 0
#define IMAGPART 1

const int MaxTrajectoryLength = 65536;
const int OutLength = MaxTrajectoryLength / 2 + 1; 

class Fourier
{
	public:
		double *inx;
		double *iny;
		double *inz;

		fftw_complex *outx;
		fftw_complex *outy;
		fftw_complex *outz;

		fftw_plan px;
		fftw_plan py;
		fftw_plan pz;

		void zero_out_input( void );
		void do_fourier( void );
		
		Fourier( void );
		~Fourier();
};

std::vector<double> fft_one_side( std::vector<double> signal );
std::vector<double> fft_two_side( std::vector<double> &signal );

void fftfreq( std::vector<double> &freqs, const int len, const double d );
void fftshift( std::vector<double> &arr, const int len );
void ifftshift( std::vector<double> &arr, const int len );

void fft_positive( std::vector<double> &signal );

