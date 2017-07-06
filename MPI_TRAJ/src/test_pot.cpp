#include "psp_pes.h"

#include <stdio.h>
#include <cstdlib>




int main()
{
	double R= 4.1;
	double Theta = 1.32;

	FILE* pfile;
	pfile = fopen("data.txt","w");

	while (R < 31)
	{
		//printf("R  = %.3f   %.12f\n",R,ab_initio_pot(R,Theta));
		fprintf(pfile,"%.2f   %.13f\n",R,dpsp_pesdTheta(R,Theta));
		R += 0.2;
 	}
 	fclose(pfile);
	system("pause");
	return 0;
}