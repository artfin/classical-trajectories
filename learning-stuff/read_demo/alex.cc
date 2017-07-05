#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main()
{
	FILE* inputfile = fopen("data.txt", "r");

	double R, Theta;

	while ( true )
	{
		int scanfResult = fwscanf(inputfile, L"%lf %lf\n", &R, &Theta);	

		if ( scanfResult == -1)
	   	{
			break;
		}

		cout << "read R: " << R << " Theta: " << Theta << endl;
	}

	fclose(inputfile);

	return 0;
}
