
//#include "const_values.h"
/* ------------------------ MODULE mgear.cpp ------------------------ */
#include "../jean-pierre/basis.h"         /*  for umleiten, fprintf, stderr, scanf,  */
                           /*       printf, NULL, REAL, LZS, LZP,     */
                           /*       fehler_melden, fehler_t           */
#include "../jean-pierre/vmblock.h"       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include "../jean-pierre/gear.h"          /*  for  gear4, gear_fehlertext            */
#include "../jean-pierre/t_dgls.h"        /*  for  bsptyp, dgls_waehlen              */

/* ------------------------------------------------------------------ */


//#include "sys/psp_pes_tapenade_d.h"



//#include "f.h"
#include "matrix.h"



/*static void syst(REAL x, REAL *y, REAL *f)
{
  //f[0] = y[0] * y[1] + COS(x) - HALF * SIN(TWO * x);
 // f[1] = y[0] * y[0] + y[1] * y[1] - (ONE + SIN(x));
  f[0] = x;
}*/


void 
syst (REAL t ,REAL *y, REAL *f)
{
  (void)(t); /* avoid unused parameter warning */
  //double * par = (double *)params;

  /*double *par = new double[4];
   par[0]= 14579;
  par[1]=38183;
  par[2]=4.398;
  par[3]=3.0;*/

  double * out = new double[6];
  rhs(out, y[0],y[1],y[2],y[3],y[4],y[5],3.0);


  f[0] = out[0];//dR/dt  
 // f[1]=f1(y,par)-(dpotdR);//d(pR)/dt 
  f[1]=out[1];
 
  f[2] = out[2];//d(Theta)/dt 
  f[3] = out[3];//d(pT)/dt
 
  f[4] = out[4];//d(phi)/dt
  f[5] = out[5];//d(theta)/dt



  delete [] out;
  //delete [] par;

  
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

/* -------------------- read input  -------------------- */

  N = 6;
 
  vmblock = vminit();                 /* initialize storage */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, N, 0);
 
  if (! vmcomplete(vmblock))  {       /* out of memory? */
    printf("mgear: out of memory.\n");
    return 0;
  }

  

  epsabs = 1E-10;
  epsrel = 1E-10;
  x0 = 0.0;

   //double y[6] = { 10.0, -1.0,1.0,1.0,1.0,1.0 };

  y0[0] = 7.0;
  y0[1] = 1.0;
  y0[2] = -0.2;
  y0[3] = 1.0;
  y0[4] = 1.0;
  y0[5] = 1.0;


  h =0.1;
  xend = 30000;
  fmax = 1000000;  /* ------------ put out the input data ----------- */



 


  /* ------------ Solve system of DEs -------------- */

  fehler = gear4(&x0, xend, N, syst, y0, epsabs,
                 epsrel, &h, fmax, &aufrufe);

  if (fehler != 0) {
    printf(" Gear4: error n° %d\n", 10 + fehler);
    return 0;
  }


  /* -------------------- put out results ------------------- */

  printf("\n\n"
         "Output data:\n"
         "------------\n"
         "error code from gear4():                  %24d\n"
         "final local step size:                    %24.15"LZP"e\n"
         "number of calls of right hand side:       %24ld\n"
         "Integration stopped at x =                %24.15"LZP"e\n\n",
         fehler, h, aufrufe, x0);

  for (i = 0; i < 1; i++)
    printf("approximate solution y%d(x) = %24.15"LZP"e\n",
           i + 1, y0[i]);

  

  return 0;
}

/* -------------------------- END mgear.cpp ------------------------- */
