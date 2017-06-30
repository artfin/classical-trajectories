#include <cmath>
#include <vector>
#include <iostream>
#include <string>

using namespace std;

// Computes te DFT of the given complex vector
void computeDFT1(const vector<double> &inreal, const vector<double> &inimag, vector<double> &outreal, vector<double> &outimag) {
	
	size_t n = inreal.size();
	for (size_t k = 0; k < n; k++) { // for each output element
		
		double sumreal = 0;
		double sumimag = 0;

		for (size_t t = 0; t < n; t++) { // for each input element
			double angle = - 2 * M_PI * t * k / n;
	
			sumreal += inreal[t] * cos(angle) + inimag[t] * sin(angle);
			sumimag += -inreal[t] * sin(angle) + inimag[t] * cos(angle);
		}

		outreal[k] = sumreal;
		outimag[k] = sumimag;
	}
}

void computeDFT2(const vector<double> &inreal, const vector<double> &inimag, vector<double> &outreal, vector<double> &outimag) {
	
	size_t n = inreal.size();
	for (size_t k = 0; k < n; k++) { // for each output element
		
		double sumreal = 0;
		double sumimag = 0;

		for (size_t t = 0; t < n; t++) { // for each input element
			double angle = - t * k / n;
	
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

void power_spectrum(const vector<double> &real, const vector<double> &imag, string name)
{
	size_t n = real.size();
	size_t power;

	cout << name << endl;

	for ( size_t i = 0; i < n; i++ ) {
		power = real[i] * real[i] + imag[i] * imag[i];

		cout << "[" << i << "]:" << power << endl;
	}

	cout << endl;
}

int main()
{
	size_t n = 8;
	static const double in1[] = {0.0, 0.707, 1.0, 0.707, 0.0, -0.707, -1.0, -0.707};
	const vector<double> in_real (in1, in1 + n);
	print_vector(in_real, "in_real");

	static const double in2[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	const vector<double> in_imag (in2, in2 + n);
	print_vector(in_imag, "in_imag");

	vector<double> out_real1(n);
	vector<double> out_imag1(n);

	vector<double> out_real2(n);
	vector<double> out_imag2(n);

	computeDFT1(in_real, in_imag, out_real1, out_imag1);
	computeDFT2(in_real, in_imag, out_real2, out_imag2);

	power_spectrum(out_real1, out_imag1, "1");
	power_spectrum(out_real2, out_imag2, "2");

	print_vector(out_real2, "out_real2");
	print_vector(out_imag2, "out_imag2");

	return 0;
}
