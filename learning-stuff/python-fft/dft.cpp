#include <cmath>
#include <vector>
#include <iostream>
#include <string>

using namespace std;

// Computes te DFT of the given complex vector
void computeDFT(const vector<double> &inreal, const vector<double> &inimag, vector<double> &outreal, vector<double> &outimag) {
	
	size_t n = inreal.size();
	for (size_t k = 0; k < n; k++) { // for each output element
		
		double sumreal = 0;
		double sumimag = 0;

		for (size_t t = 0; t < n; t++) { // for each input element
			double angle = 2 * M_PI * t * k / n;
	
			sumreal += inreal[t] * cos(angle) + inimag[t] * sin(angle);
			sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
		}

		outreal[k] = sumreal;
		outimag[k] = sumimag;
	}
}

void print_vector(const vector<double> &v, string name) {

	size_t n = v.size();

	for (size_t i = 0; i < n; i++) {
		cout << name << "[" << i << "]: " << v[i] << endl;
	}
	cout << "=============" << endl;
}

int main()
{
	size_t n = 4;
	static const double arr1[] = {1.0, 2.0, 3.0, 4.0};
	const vector<double> in_real (arr1, arr1 + n);
	print_vector(in_real, "in_real");

	static const double arr2[] = {0.0, 0.0, 0.0, 0.0};
	const vector<double> in_imag (arr2, arr2 + n);
	print_vector(in_imag, "in_imag");

	vector<double> out_real(4);
	vector<double> out_imag(4);

	computeDFT(in_real, in_imag, out_real, out_imag);

	print_vector(out_real, "out_real");
	print_vector(out_imag, "out_imag");

	return 0;
}
