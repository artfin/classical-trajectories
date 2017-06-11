/* ------------------------ MODULE mgear.cpp ------------------------ */
#include "../jean-pierre/basis.h"         /*  for umleiten, fprintf, stderr, scanf,  */
                           		  /*       printf, NULL, REAL, LZS, LZP,     */
                           		  /*       fehler_melden, fehler_t           */
#include "../jean-pierre/vmblock.h"       /*  for  vmalloc, vmcomplete, vmfree,      */
                          		  /*       vminit, VEKTOR                    */
#include "../jean-pierre/gear.h"          /*  for  gear4, gear_fehlertext            */
#include "../jean-pierre/t_dgls.h"        /*  for  bsptyp, dgls_waehlen              */
/* ------------------------------------------------------------------ */

#include "matrix.h"
#include <iostream>

using namespace std;

const double J = 5.0;

void syst (REAL t ,REAL *y, REAL *f)
{
  (void) (t);  /* avoid unused parameter warning */

  double * out = new double[6];
  rhs(out, y[0], y[1], y[2], y[3], y[4], y[5], J);

  f[0] = out[0];   // dR/dt  
  f[1] = out[1];   // d(pR)/dt
  f[2] = out[2];   // d(Theta)/dt 
  f[3] = out[3];   // d(pT)/dt
  f[4] = out[4];   // d(phi)/dt
  f[5] = out[5];   // d(theta)/dt

  delete [] out;
}


int main() {
 
  REAL     epsabs;       /* absolute error bound                      */
  REAL     epsrel;       /* relative error bound                      */
  REAL     x0;           /* left edge of integration interval         */
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

  N = 6;

  // initialize storage
  vmblock = vminit();           
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, N, 0);
  
  // out of memory?
  if (! vmcomplete(vmblock))  {
    printf("mgear: out of memory.\n");
    return 0;
  }

  epsabs = 1E-10;
  epsrel = 1E-10;
  
  int T_STEP = 100;
  int T_LENGTH= 5500000;

  h = 0.1;
  fmax = 100000;  

  // y = [R, pR, theta, pT, alpha, beta]
  y0[0] = 20.0;
  y0[1] = -5.0;
  y0[2] = -0.2;
  y0[3] = 1.0;
  y0[4] = 1.0;
  y0[5] = 1.0;

  FILE *traj_file = fopen("output/trajectory.dat", "w");
  FILE *dipole_file = fopen("output/dipole.dat", "w");

  for ( int T = 0; T < T_LENGTH - 1; T += T_STEP ) {
      
	x0 = T;
      	xend = T + T_STEP;
	cout << "x0: " << x0 << "; xend: " << xend << endl;
 
  	fehler = gear4(&x0, xend, N, syst, y0, epsabs, epsrel, &h, fmax, &aufrufe);

  	double *dip_out = new double[4];
  	hamiltonian(dip_out, y0[0], y0[1], y0[2], y0[3], y0[4], y0[5], J, true);

	if ( dipole_file != NULL ) {
		fprintf(dipole_file, "%lf %6lf %6lf %6lf %6lf\n", xend, dip_out[0], dip_out[1], dip_out[2], dip_out[3]);
	}

	if ( traj_file != NULL ) {
		fprintf(traj_file, "%lf %6lf %6lf %6lf %6lf %6lf %6lf %6lf\n", xend, y0[0], y0[1], y0[2], y0[3], y0[4], y0[5], J);
	}
	
  	if (fehler != 0) {
            printf(" Gear4: error n° %d\n", 10 + fehler);
            return 0;
  	}
  } 
  
  fclose(traj_file);
  fclose(dipole_file);

  return 0;
}

