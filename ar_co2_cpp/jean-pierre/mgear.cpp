/* ------------------------ MODULE mgear.cpp ------------------------ */
#include "basis.h"         /*  for umleiten, fprintf, stderr, scanf,  */
                           /*       printf, NULL, REAL, LZS, LZP,     */
                           /*       fehler_melden, fehler_t           */
#include "vmblock.h"       /*  for  vmalloc, vmcomplete, vmfree,      */
                           /*       vminit, VEKTOR                    */
#include "gear.h"          /*  for  gear4, gear_fehlertext            */
#include "t_dgls.h"        /*  for  bsptyp, dgls_waehlen              */

/* ------------------------------------------------------------------ */

int main() {
/***********************************************************************
*    Solve a first order Stiff System of Differential Equations using  *
*    the implicit Gear method of order 4.                              *
* -------------------------------------------------------------------  *
* Mode of operation:                                                   *
* ==================                                                   *
* This program solves one of the 11 examples of file t_dgls.cpp using  *
* the implicit Gear method of order 4 (see file gear.cpp).             *
* To test other systems of DEs, please proceed as explained in file    *
* t_dgls.cpp.                                                          *
*                                                                      *
*   Inputs:                                                            *
*   =======                                                            *
*   bspnummer  Number of DE system from t_dgls.cpp                     *
*   epsabs     desired absolute error bound                            *
*   epsrel     desired relative error bound                            *
*   x0         left edge of integration                                *
*   y0[0]   \  known approximation for the solution at x0              *
* ..  .      >                                                         *
*   y0[n-1] /                                                          *
*   h          initial step size                                       *
*   xend       right endpoint of integration                           *
*   fmax       maximal number of calls of the right hand side          *
*                                                                      *
*   The size n of the DE system is passed on from t_dgls.cpp.          *
* -------------------------------------------------------------------- *
* SAMPLE RUN                                                           *
*                                                                      *
* Example #1:                                                          *
* (Solve set of differential equations (n=2):                          *
*     f[0] = y[0] * y[1] + COS(x) - HALF * SIN(TWO * x);               *
*     f[1] = y[0] * y[0] + y[1] * y[1] - (ONE + SIN(x));               *
*  Find values of f(0), f(1) at x=1.5).                                *
*                                                                      *
* Input example number (0 to 11): 0                                    *
* abs. epsilon: 1e-6                                                   *
* rel. epsilon: 1e-8                                                   *
* x0: 0                                                                *
* y0[0]: 0.5                                                           *
* y0[1]: 0.5                                                           *
* initial step size h: 0.0001                                          *
* right edge xend: 1.5                                                 *
* maximal number of calls of right hand side: 6000                     *
*                                                                      *
* Input data:                                                          *
* -----------                                                          *
* Example  =                       0                                   *
* n        =                       2                                   *
* x0       =  0.000000000000000e+000                                   *
* xend     =  1.500000000000000e+000                                   *
* epsabs   =  1.000000000000000e-006                                   *
* epsrel   =  1.000000000000000e-008                                   *
* fmax     =                    6000                                   *
* h        =  1.000000000000000e-004                                   *
* y0[0]    =  5.000000000000000e-001                                   *
* y0[1]    =  5.000000000000000e-001                                   *
*                                                                      *
* Output data:                                                         *
* ------------                                                         *
* error code from gear4():                                 0           *
* final local step size:              6.067837631351808e-002           *
* number of calls of right hand size:                    355           *
* Integration stopped at x =          1.500000000000000e+000           *
*                                                                      *
* approximate solution y1(x) =   1.235986128938243e+000                *
* approximate solution y2(x) =  -1.049496179812494e-001                *
*                                                                      *
* Example #2:                                                          *
* (Solve set of differential equations (n=5):                          *
*   f[0] = y[1];                                                       *
*   f[1] = y[2];                                                       *
*   f[2] = y[3];                                                       *
*   f[3] = y[4];                                                       *
*   f[4] = ((REAL)45.0 * y[2] * y[3] * y[4] -                          *
          (REAL)40.0 * y[3] * y[3] * y[3]) / (NINE * y[2] * y[2]);     *
*  Find values of f(0), ..., f(4) at x=1.5).                           *
*                                                                      *
* Input example number (0 to 11): 3                                    *
* abs. epsilon: 1e-10                                                  *
* rel. epsilon: 1e-10                                                  *
* x0: 0                                                                *
* y0[0]: 1                                                             *
* y0[1]: 1                                                             *
* y0[2]: 1                                                             *
* y0[3]: 1                                                             *
* y0[4]: 1                                                             *
* initial step size h: 0.001                                           *
* right edge xend: 1.5                                                 *
* maximal number of calls of right hand side: 6000                     *
*                                                                      *
* Input data:                                                          *
* -----------                                                          *
* Example  =                       3                                   *
* n        =                       5                                   *
* x0       =  0.000000000000000e+000                                   *
* xend     =  1.500000000000000e+000                                   *
* epsabs   =  1.000000000000000e-010                                   *
* epsrel   =  1.000000000000000e-010                                   *
* fmax     =                    6000                                   *
* h        =  1.000000000000000e-003                                   *
* y0[0]    =  1.000000000000000e+000                                   *
* y0[1]    =  1.000000000000000e+000                                   *
* y0[2]    =  1.000000000000000e+000                                   *
* y0[3]    =  1.000000000000000e+000                                   *
* y0[4]    =  1.000000000000000e+000                                   *
*                                                                      *
* Output data:                                                         *
* ------------                                                         *
* error code from gear4():                                 0           *
* final local step size:              4.863476630190404e-003           *
* number of calls of right hand size:                   3418           *
* Integration stopped at x =          1.500000000000000e+000           *
*                                                                      *
* approximate solution y1(x) =   4.363961029902760e+000                *
* approximate solution y2(x) =   4.000000007634299e+000                *
* approximate solution y3(x) =   2.828427156619908e+000                *
* approximate solution y4(x) =   4.861628881847740e-008                *
* approximate solution y5(x) =  -3.771236222955674e+000                *
*                                                                      *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges        *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
*                                                                      *
*                              C++ Release By J-P Moreau, Paris.       *
*                                     (www.jpmoreau.fr)                *
***********************************************************************/
// Add to project the following files:
// ==================================
// awp.h, awp.cpp, basis.h, basis_r.cpp, fgauss.cpp, gear.h, gear.cpp,
// t_dgls.h,t_dgls.cpp, vmbloch.h, vmblock.cpp.

// NOTE: Utility files basis.h, basis_r.cpp, vmblock.h, vmblock.cpp are
// located in directory MATRICES. Other files are in directory DIFFEQUA.
 
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

  printf("\n Input example number (0 to 11): ");
  scanf("%d", &bspnummer);
  if ((beispiel = dgls_waehlen(bspnummer)) == NULL) {
    printf("\n tdgls: non registered example.\n");
    return 0;
  }

  n = beispiel->n;
  vmblock = vminit();                 /* initialize storage */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  yex = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))  {       /* out of memory? */
    printf("mgear: out of memory.\n");
    return 0;
  }

  printf(" abs. epsilon: ");
  scanf("%"LZS"f", &epsabs);

  printf(" rel. epsilon: ");
  scanf("%"LZS"f", &epsrel);

  printf(" x0: ");
  scanf("%"LZS"f", &x0);

  for (i = 0; i < n; i++)
  {
    printf(" y0[%d]: ",i);
    scanf("%"LZS"f", y0 + i);
  }

  printf(" initial step size h: ");
  scanf("%"LZS"f", &h);

  printf(" right edge xend: ");
  scanf("%"LZS"f", &xend);

  printf(" maximal number of calls of right hand side: ");
  scanf("%ld", &fmax);


  /* ------------ put out the input data ----------- */

  printf("\n"
         "Solve a first order ordinary system of DEs\n"
         "==========================================\n"
         "using the implicit method of Gear of 4th order\n"
         "==============================================\n\n\n"
         "System of DEs:\n"
         "--------------\n"
         "%s\n\n"
         "Input data:\n"
         "-----------\n"
         "Example  = %24d\n"
         "n        = %24d\n"
         "x0       = %24.15"LZP"e\n"
         "xend     = %24.15"LZP"e\n"
         "epsabs   = %24.15"LZP"e\n"
         "epsrel   = %24.15"LZP"e\n"
         "fmax     = %24ld\n"
         "h        = %24.15"LZP"e\n",
         (*beispiel->dgl_text)(), bspnummer, n, x0, xend, epsabs,
         epsrel, fmax, h);

  for (i = 0; i < n; i++)
    printf("y0[%d]    = %24.15"LZP"e\n", i, y0[i]);


  /* ------------ Solve system of DEs -------------- */

  fehler = gear4(&x0, xend, n, beispiel->rechte_seite, y0, epsabs,
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

  for (i = 0; i < n; i++)
    printf("approximate solution y%d(x) = %24.15"LZP"e\n",
           i + 1, y0[i]);

  if (beispiel->exakte_loesung != NULL)       /* "exact" solution     */
  {                                           /* available?           */
    (*beispiel->exakte_loesung)(x0, yex);
    printf("\n");
    for (i = 0; i < n; i++)
      printf("'exact' solution     y%d(x) = %24.15"LZP"e\n",
             i + 1, yex[i]);
    printf("\nDifference  approximate solution - 'exact' solution:\n");
    for (i = 0; i < n; i++)
      printf("%24.15"LZP"g\n", y0[i] - yex[i]);
  }

  return 0;
}

/* -------------------------- END mgear.cpp ------------------------- */
