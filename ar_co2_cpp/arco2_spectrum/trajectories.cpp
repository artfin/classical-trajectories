/* ------------------------ MODULE mgear.cpp ------------------------ */
#include "../jean-pierre/basis.h"         
/*  for umleiten, fprintf, stderr, scanf,  */
/*       printf, NULL, REAL, LZS, LZP,     */
/*       fehler_melden, fehler_t           */
#include "../jean-pierre/vmblock.h"       
/*  for  vmalloc, vmcomplete, vmfree,      */
/*       vminit, VEKTOR                    */
#include "../jean-pierre/gear.h"          
/*  for  gear4, gear_fehlertext            */
#include "../jean-pierre/t_dgls.h"        
/*  for  bsptyp, dgls_waehlen              */
/* ------------------------------------------------------------------ */

#include "matrix.h"
#include <vector>
#include <fftw3.h>

using namespace std;

// defining total angular momentum J
#define ANG 1

// macros for real and imaginary parts
#define REALPART 0
#define IMAGPART 1

// Computes te DFT of the given complex vector
// FFT in fftw3 computes DFT with angle = - 2 * M_PI * t * k / n
void computeDFT2(const vector<double> &inreal, vector<double> &outreal, vector<double> &outimag) {
	
	size_t n = inreal.size();
	for (size_t k = 0; k < n; k++) { // for each output element
		
		double sumreal = 0;
		double sumimag = 0;

		for (size_t t = 0; t < n; t++) { // for each input element
			double angle = - t * k / n;
	
			sumreal += inreal[t] * cos(angle);
			sumimag += -inreal[t] * sin(angle);
		}

		outreal[k] = sumreal;
		outimag[k] = sumimag;
	}
}

void syst (REAL t ,REAL *y, REAL *f)
{
  (void)(t); // avoid unused parameter warning 

  double *out = new double[6];
  rhs(out, y[0], y[1], y[2], y[3], y[4], y[5], ANG);
  // R  Theta pR pT phi theta J

  f[0] = out[0]; // dR/dt  
  f[1] = out[1]; // d(Theta)/dt
  f[2] = out[2]; // d(pR)/dt
  f[3] = out[3]; // d(pT)/dt
  f[4] = out[4]; // d(phi)/dt
  f[5] = out[5]; // d(theta)/dt

  delete [] out;
}

int main() {

  REAL     epsabs;       /* absolute error bound                      */
  REAL     epsrel;       /* relative error bound                      */
  REAL     t0;           /* left edge of integration interval         */
  REAL     *y0;          /* [0..n-1]-vector: initial value, approxim. */
        
  REAL     h;            /* initial, final step size                  */
  REAL     xend;         /* right edge of integration interval        */
  long     fmax;         /* maximal number of calls of right side     */
                         /* in gear4()                                */
  long     aufrufe;      /* actual number of function calls           */

  int      N;            /* number of DEs in system                   */
  int      fehler;       /* error code from umleiten(), gear4()       */
  int      i;            /* loop counter                              */
                         
  void     *vmblock;     /* List of dynamically allocated vectors     */

  FILE *trajectory_file = fopen("_traj.txt", "w");
  FILE *dipfft = fopen("dipfft.txt", "w"); 
  FILE *testfft1 = fopen("testfft1.txt", "w");
  FILE *testfft2 = fopen("testfft2.txt", "w");
  /* -------------------- read input  -------------------- */

  N = 6;
 
  vmblock = vminit();                 /* initialize storage */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, N, 0);
 
  if (! vmcomplete(vmblock))  {       /* out of memory? */
    printf("mgear: out of memory.\n");
    return 0;
  }

  double step = 3000;
  double end = 700000;

  int num_points = (end / step);

  epsabs = 1E-13;
  epsrel = 1E-13;
  
  t0 = 0.0;
 
  y0[0] = 15.0;
  y0[1] = 1.0;
  y0[2] = -4.0;
  y0[3] = 1.0;
  y0[4] = 1.0;
  y0[5] = 1.0;

  h = 0.1;         // initial step size
  xend = step;     // initial right bound of integration
  fmax = 1000000;  // maximal number of calls 

  double *dipole = new double [3];

  vector<double> ddipx(num_points);
  vector<double> ddipy(num_points);
  vector<double> ddipz(num_points);

  for (int i = 0; i < num_points; ++i)
  {
     fehler = gear4(&t0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);
     
     if ( fehler != 0 ) {
     	printf(" Gear4: error n° %d\n", 10 + fehler);
     	return 0;
     }
     
     hamiltonian(dipole, y0[0], y0[1], y0[2], y0[3], y0[4], y0[5], ANG, true);
     
     //fprintf(trajectory_file, "%f %.12f %d\n", t0, y0[0], aufrufe);
     fprintf(trajectory_file, "%f %.12f %.12f %.12f %.12f %.12f\n", t0, y0[0], y0[1], dipole[0], dipole[1], dipole[2]);

     ddipx.push_back(dipole[0]);
     ddipy.push_back(dipole[1]);
     ddipz.push_back(dipole[2]);

     xend = step * (i + 2);
     aufrufe = 0;  // actual number of calls
  }
 
  fclose(trajectory_file);

  // length of dipole vector
  size_t n = ddipx.size();

  // input and output arrays
  fftw_complex _ddipx[n];
  fftw_complex _ddipx_fftw[n];

  fftw_complex _ddipy[n];
  fftw_complex _ddipy_fftw[n];

  fftw_complex _ddipz[n];
  fftw_complex _ddipz_fftw[n];

  // filling arrays for fftw
  for ( int i = 0; i < n; i++ )
  {
	  _ddipx[i][REALPART] = ddipx[i];
	  _ddipx[i][IMAGPART] = 0; 

	  _ddipy[i][REALPART] = ddipy[i];
	  _ddipy[i][IMAGPART] = 0; 

	  _ddipz[i][REALPART] = ddipz[i];
	  _ddipz[i][IMAGPART] = 0;
  }

  // planning the FFT and executing it
  fftw_plan planx = fftw_plan_dft_1d(n, _ddipx, _ddipx_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan plany = fftw_plan_dft_1d(n, _ddipy, _ddipy_fftw, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan planz = fftw_plan_dft_1d(n, _ddipz, _ddipz_fftw, FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(planx);
  fftw_execute(plany);
  fftw_execute(planz);

  // do some cleaning
  fftw_destroy_plan(planx);
  fftw_destroy_plan(plany);
  fftw_destroy_plan(planz);
  
  fftw_cleanup(); 

  //vector<double> _ddipx_fft2_real(n);
  //vector<double> _ddipx_fft2_imag(n);
  //vector<double> _ddipy_fft2_real(n);
  //vector<double> _ddipy_fft2_imag(n);
  //vector<double> _ddipz_fft2_real(n);
  //vector<double> _ddipz_fft2_imag(n);

  //computeDFT2(ddipx, _ddipx_fft2_real, _ddipx_fft2_imag);
  //computeDFT2(ddipy, _ddipy_fft2_real, _ddipy_fft2_imag);
  //computeDFT2(ddipz, _ddipz_fft2_real, _ddipz_fft2_imag);

  // px -- square amplitude (power) of x-component
  // py -- square amplitude (power) of y-component
  // pz -- square amplitude (power) of z-component
  // p  -- sum of square amplitudes of all components  
  double p, px, py, pz;
 
  double p2, px2, py2, pz2;
  
  for ( int i = 0; i < n; i++ )
  {
	  px = _ddipx_fftw[i][REALPART] * _ddipx_fftw[i][REALPART] + _ddipx_fftw[i][IMAGPART] * _ddipx_fftw[i][IMAGPART];
	  py = _ddipy_fftw[i][REALPART] * _ddipy_fftw[i][REALPART] + _ddipy_fftw[i][IMAGPART] * _ddipy_fftw[i][IMAGPART];
	  pz = _ddipz_fftw[i][REALPART] * _ddipz_fftw[i][REALPART] + _ddipz_fftw[i][IMAGPART] * _ddipz_fftw[i][IMAGPART];

	  //px2 = _ddipx_fft2_real[i] * _ddipx_fft2_real[i] + _ddipx_fft2_imag[i] * _ddipx_fft2_imag[i];
	  //py2 = _ddipy_fft2_real[i] * _ddipy_fft2_real[i] + _ddipy_fft2_imag[i] * _ddipy_fft2_imag[i];
	  //pz2 = _ddipz_fft2_real[i] * _ddipz_fft2_real[i] + _ddipz_fft2_imag[i] * _ddipz_fft2_imag[i];	  
	  
	  p = px + py + pz;
	  //p2 = px2 + py2 + pz2;

	  fprintf(dipfft, "%.12f\n", p);
	  //fprintf(dipfft, "%.12f %.12f\n", p, p2);
  }

  fclose(dipfft);

  return 0;
}

/* -------------------------- END mgear.cpp ------------------------- */
