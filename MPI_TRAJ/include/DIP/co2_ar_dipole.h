#pragma once 

#include <cmath>

double dipx(double R, double Theta);
double dipz(double R, double Theta);

double dipx_cos(double R_inv,double cos_theta);
double dipz_cos(double R_inv,double cos_theta);

double ddipxdR(double R, double Theta);
double ddipzdR(double R, double Theta);

double ddipxdTheta(double R, double Theta);
double ddipzdTheta(double R, double Theta);
