/* ------------------------ MODULE mgear.cpp ------------------------ */
#include "basis.h"         /*  for umleiten, fprintf, stderr, scanf,  */
                           /*       printf, NULL, REAL, LZS, LZP,     */
                           /*       fehler_melden, fehler_t           */
#include "vmblock.h"       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include "gear.h"          /*  for  gear4, gear_fehlertext            */
#include "t_dgls.h"        /*  for  bsptyp, dgls_waehlen              */

/* ------------------------------------------------------------------ */


#include "sys/psp_pes_tapenade_d.h"

#include "sys/f.h"



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

  double *par = new double[4];
   par[0]= 14579;
  par[1]=38183;
  par[2]=4.398;
  par[3]=3.0;

  double  dpotdR;
  double  dpotdTheta;
  double result;
  potential_d(y[0],1.0,y[2],0.0,result,dpotdR);
  potential_d(y[0],0.0,y[2],1.0,result,dpotdTheta);

  f[0] = f4(y,par);//dR/dt  
  f[1]=f1(y,par)-(dpotdR);//d(pR)/dt 
 
  f[2] = f3(y,par);//d(Theta)/dt 
  f[3] = f2(y,par)-(dpotdTheta);//d(pT)/dt
 
  f[4] = f5(y,par);//d(phi)/dt
  f[5] = f6(y,par);//d(theta)/dt
  delete [] par;

  
}




int main() {

 
  REAL     epsabs;       /* absolute error bound                      */
  REAL     epsrel;       /* relative error bound                      */
  REAL     x0;           /* left edge of integration interval         */
  REAL     *y0;          /* [0..n-1]-vector: initial value, approxim. */
  REAL     *yex;         /* [0..n-1]-vector: exact solution           */
  REAL     h;            /* initial, final step size                  */
  REAL     xend;         /* right edge of integration interval        */
  long     fmax;         /* maximal number of calls of right side     */
                         /* in gear4()                                */
  long     aufrufe;      /* actual number of function calls           */
  int      bspnummer;    /* Number of the system of DEs from t_dgls.c */
  int      n;            /* number of DEs in system                   */
  int      fehler;       /* error code from umleiten(), gear4()       */
  int      i;            /* loop counter                              */
  bsptyp   *beispiel;    /* pointer to the structure that describes   */
                         /* the actual system of DEs                  */
  void     *vmblock;     /* List of dynamically allocated vectors     */

/* -------------------- read input  -------------------- */

  n = 6;
 
  vmblock = vminit();                 /* initialize storage */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  yex = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))  {       /* out of memory? */
    printf("mgear: out of memory.\n");
    return 0;
  }

  

  epsabs = 1E-9;
  epsrel = 1E-9;
  x0 = 0.0;

   //double y[6] = { 10.0, -1.0,1.0,1.0,1.0,1.0 };

  y0[0] = 7.0;
  y0[1] = -0.2;
  y0[2] = 1.0;
  y0[3] = 1.0;
  y0[4] = 1.0;
  y0[5] = 1.0;


  h =0.1;
  xend = 30000;
  fmax = 1000000;  /* ------------ put out the input data ----------- */



 


  /* ------------ Solve system of DEs -------------- */

  fehler = gear4(&x0, xend, n, syst, y0, epsabs,
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
