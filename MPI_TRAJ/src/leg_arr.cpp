#include "leg_arr.h"

double * legendre_array(int N, double cosT)
{
	double * leg = (double*)malloc((N+1)*sizeof(double));
	leg[0] = 1;	
	leg[1] = cosT;

	for (int i =1; i < N; ++i)
	{
		leg[i+1] = ((2.0*i+1.0)/(i+1.0))*cosT*leg[i] - ((double)i/(i+1.0))*leg[i-1];
	}
	return leg;

}
