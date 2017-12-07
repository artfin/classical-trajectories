/* ------------------------ MODULE gear.cpp ------------------------- */

/***********************************************************************
*                                                                      *
* Solve a first order system of DEs using the implicit Gear method     *
* of order 4.                                                          *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Klaus Niederdrenk, 1.22.1996 (FORTRAN 77)      *
* Adaptation:           Juergen Dietel, Computing Center, RWTH Aachen  *
* Source:               FORTRAN 77 source code                         *
* Date:                 2.26.1996                                      *
*                                                                      *
*                       C++ Release By J-P Moreau, Paris.              *
*                       (www.jpmoreau.fr)                              *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges        *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
***********************************************************************/

#include "basis.h"     /*  for  REAL, ZERO, FABS, dglsysfnk, POW,     */
                       /*       MACH_EPS, TEN, SQRT, ONE, norm_max,   */
                       /*       copy_vector, FIVE, TWO, THREE, FOUR,  */
                       /*       gauss, NULL, SIX, sqr, min,           */
                       /*       max, boolean, TRUE, FALSE, fehler_t   */
#include "vmblock.h"   /*  for  vmalloc, vmcomplete, vmfree, vminit,  */
                       /*       VEKTOR, MATRIX, VVEKTOR               */
#include "awp.h"       /*  for  awp                                   */
#include "gear.h"      /*  for  gear4, gear_fehlertext                */

/* ------------------------------------------------------------------ */

#ifdef DEBUG
/* put put a REAL variable `var' (and its name)                 */
#define zeigrvar(var)  printf("%-12s%16"LZP"g\n", #var": ", var);

/* put out a  REAL vector of length `n' (and its name)       */
#define zeigrvek(vek, n)              \
  {                                   \
    int jjj;                          \
    printf("%-12s", #vek": ");        \
    for (jjj = 0; jjj < n; jjj++)     \
      printf("%16"LZP"g", vek[jjj]);  \
    printf("\n");                     \
  }

/* put out a REAL square matrix with name                   */
#define zeigrmat(mat, n)                          \
{                                                 \
  int iii, jjj;                                   \
  for (iii = 0; iii < (n); iii++)                 \
  {                                               \
    printf("%-12s", (iii == 0) ? #mat": " : "");  \
    for (jjj = 0; jjj < (n); jjj++)               \
      printf("%16"LZP"g", (mat)[iii][jjj]);       \
    printf("\n");                                 \
  }                                               \
}

#else
#define zeigrvar(var)
#define zeigrvek(vek, n)
#define zeigrmat(mat, n)
#endif



/* ------------------------------------------------------------------ */

static REAL dist_max   /* Maximum norm of difference of two vectors...*/
/*.IX{dist\unt max}*/
    (
     REAL      vektor1[],
     REAL      vektor2[],
     int       n
    )

/***********************************************************************
* Find distance of vectors vektor1 and vektor2 in maximum norm.        *
*                                                                      *
* global variables:                                                    *
* =================                                                    *
* REAL, ZERO, FABS                                                     *
***********************************************************************/

{
  REAL abstand,         /* aux value */
       hilf;            /* ditto     */

  for (n--, abstand = ZERO; n >= 0; n--)
    if ((hilf= FABS(vektor1[n] - vektor2[n])) > abstand)
      abstand = hilf;

  return abstand;
}



/*--------------------------------------------------------------------*/

int gear4          /* Gear  method of 4th order for DESs of 1st order */
    (
     REAL      *x,             /* starting or end point      .........*/
     REAL      xend,           /* desired end point (> x)     ........*/
     int       n,              /* number of DEs    ...................*/
	 //std::function<void()> dgl,
     dglsysfnk dgl,            /* right hand side of system of DEs ...*/
     REAL      y[],            /* initial value or solution ..........*/
     REAL      epsabs,         /* absolute error bound    ............*/
     REAL      epsrel,         /* relative error bound    ............*/
     REAL      *h,             /* starting or final step size   ......*/
     long      fmax,           /* maximal number of calls of dgl() ...*/
     long      *aufrufe        /* actual number of calls of dgl() ....*/
    )                          /* error code .........................*/

/***********************************************************************
* Compute the value of the solution at xend, starting with the IVP.    *
* We use the implicit method of Gear of 4th order with step size       *
* control which is especially suited for stiff DEs.                    *
* The local step size control insures that the two error bounds are met*
* The number of function calls of the right hand side is limited by    *
* fmax. This function can be used inside a loop to find solutions at   *
* a specified point to arbitrary accuracy.                             *
.BE*)
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        initial value x0                                            *
* xend     final value for the integration (xend > x0)                 *
* n        number of DEs                                               *
* dgl      Function to compute the right hand side f(x0,y0) for the    *
*          system of DEs                                               *
* y        [0..n-1] solution vector y0 of the system of DEs at x0      *
* epsabs   absolute error bound (>= 0); if zero, we only check for the *
*          relative error.                                             *
* epsrel   relative error bound (>= 0); if zero, we only check for the *
*          absolute error.                                             *
* h        given step size for first step                              *
* fmax     maximal number of calls of the right hand side of the system*
*                                                                      *
* Output parameters:                                                   *
* =================                                                    *
* x        final x-value of the integration; normally equal to xend    *
* h        final step size; keep for further calls                     *
* y        [0..n-1]-vector, the solution of the system of DEs at x     *
* aufrufe  counter of calls of dgl()                                   *
*                                                                      *
* Return values:                                                       *
* ==============                                                       *
* Error code:                                                          *
*   =   0:  all ok                                                     *
*   =   1:  Both error bounds too small                                *
*   =   2:  xend <= x0                                                 *
*   =   3:  Step size h <= 0                                           *
*   =   4:  n <= 0                                                     *
*   =   5:  # calls > fmax; we have not reached the desired xend;      *
*           repeat function call with actual values of x, y, h.        *
*   =   6:  Jacobi matrix is singular; x, y, h are the last values     *
*           that could be computed with accuracy                       *
*   =   7:  ran out of memory                                          *
*   =   8:  error when calling gauss() for the second time;            *
*           should not occur.                                          *
*                                                                      *
* global variables:                                                    *
* =================                                                    *
* REAL, POW, MACH_EPS, TEN, SQRT, ZERO, ONE, norm_max, FABS, vminit,   *
* vmalloc, VEKTOR, MATRIX, VVEKTOR, vmcomplete, vmfree, copy_vector,   *
* FIVE, TWO, THREE, FOUR, gauss, NULL, SIX, dist_max, sqr, min, max,   *
* boolean, TRUE, FALSE, awp                                            *
.BA*)
***********************************************************************/
{
  /* ---------------------- static variables --------------------- */
  static
    REAL eps1;            /* MACH_EPS ^0.75; used instead of MACH_EPS */
                          /* to check whether we have reached xend in */
                          /* order to avoid miniscule step sizes      */
  static
    REAL eps2;            /* 100 * MACH_EPS; for comparison with zero */
  static
    REAL hs;              /* optimal step size for Jacobi matrix      */
                                                /* approximation      */
  static
    boolean erster_aufruf = TRUE;                     /* first call?  */

  /* --------------------- dynamic variables --------------------- */
  void   *vmblock;             /* List of dynamic vectors and matrices*/
  REAL   *hilf;                /* [0..n-1]-vector                     */
  REAL   **zj;                 /* [0..4,0..n-1]-matrix                */
  REAL   **zjp1;               /* [0..4,0..n-1]-matrix                */
  REAL   *f;                   /* [0..n-1]-vector                     */
  REAL   *ykp1;                /* [0..n-1]-vector                     */
  REAL   **fs;                 /* [0..n-1,0..n-1]-matrix              */
  REAL   **fsg;                /* [0..n-1,0..n-1]-matrix              */
  REAL   *con;                 /* [0..n-1]-vector                     */
  int    *perm;                /* [0..n-1] permutation vector for     */
                               /* Gauss elimination                   */

  /* -------------------- automatic variables -------------------- */
  int     fehler;              /* error code                          */
  REAL    sg;                  /* sign of xend                        */
  REAL    xe;                  /* |xend| - eps2, carrying the sign    */
                               /* of xend                             */
  int     amEnde;              /* Flag, indicating that we shall reach*/
                               /* xend with the current step          */
  REAL    ymax;                /* Maximum norm of the current         */
                               /* approximate value of y              */
  REAL    dummy = ZERO;        /* aux storage for h                   */
  REAL    xka;
  REAL    xke;
  REAL    hka;
  REAL    hk1;
  REAL    diff;
  REAL    eps;
  REAL    q;
  REAL    halt;
  REAL    quot1;
  REAL    quot2;
  REAL    quot3;
  REAL    quot4;
  boolean nochmal;
  long    aufrufe_awp;
  int     signdet;             /* sign of determinant in Gauss        */
  int     iter;
  int     i;
  int     k;


  /* ------------------------ Initialize ------------------------ */

  if (erster_aufruf)
  {
    eps1          = POW(MACH_EPS, (REAL)0.75);
    eps2          = (REAL)100 * MACH_EPS;
    hs            = TEN * SQRT(MACH_EPS);
    erster_aufruf = FALSE;
  }
  zeigrvar(MACH_EPS); zeigrvar(eps1); zeigrvar(eps2); zeigrvar(hs);

  sg       = (xend >= ZERO) ? ONE : -ONE;
  xe       = (ONE - sg * eps2) * xend;
  fehler   = 0;
  *aufrufe = 0l;
  amEnde   = 0;

  /* -------- check input parameters ------- */

  ymax = norm_max(y, n);
  if (epsabs <= eps2 * ymax && epsrel <= eps2)
    fehler = 1;
  else if (xe < *x)
    fehler = 2;
  else if (*h < eps2 * FABS(*x))
    fehler = 3;
  else if (n <= 0)
    fehler = 4;
  if (fehler)
    return fehler;

  /* ------------ allocate dynamic vectors and matrices ---------- */

  vmblock = vminit();
  hilf = (REAL *) vmalloc(vmblock, VEKTOR,  n+1, 0);
  f    = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  ykp1 = (REAL *) vmalloc(vmblock, VEKTOR,  n+1, 0);
  con  = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);
  perm = (int  *) vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));
  zj   = (REAL **)vmalloc(vmblock, MATRIX,  5, n);
  zjp1 = (REAL **)vmalloc(vmblock, MATRIX,  5, n);
  fs   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  fsg  = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  if (! vmcomplete(vmblock))     /* lack of memory?    */
  {
    vmfree(vmblock);             /* free allocated storage            */
    return 7;                    /* return error                      */
  }

  /* ----------- first integration step ---------- */

  if (*x + *h > xe)
  {
    *h     = xend - *x;
    dummy  = *h;
    amEnde = 1;
  }
  copy_vector(hilf, y, n+1);
  xka = *x;
  xke = xka;
  hka = (REAL)0.25 * *h;
  hk1 = hka;
  for (k = 1; k < 5; k++)
  {
    xke += hka;
    fehler = awp(&xka, xke, n, dgl, hilf, epsabs, epsrel, &hk1,
                 6, fmax - *aufrufe, &aufrufe_awp);
    zeigrvek(hilf, n);
    *aufrufe += aufrufe_awp;
    if (fehler)
    {
      vmfree(vmblock);
      return fehler;
    }
    copy_vector(zjp1[k], hilf, n);
  }
  dgl(*x, y, f);
  ++(*aufrufe);

  /* --------- Compute first Gear-Nordsieck approximation -----*/

  for (i = 0; i < n; i++)
  {
    zj[0][i] = y[i];
    zj[1][i] = *h * f[i];
    zj[2][i] = ONE / (REAL)24.0 * ((REAL)35.0 * y[i] -
               (REAL)104.0 * zjp1[1][i] +
               (REAL)114.0 * zjp1[2][i] -
               (REAL)56.0  * zjp1[3][i] +
               (REAL)11.0  * zjp1[4][i]);
    zj[3][i] = ONE / (REAL)12.0 * (-FIVE * y[i] +
               (REAL)18.0  * zjp1[1][i] -
               (REAL)24.0  * zjp1[2][i] +
               (REAL)14.0  * zjp1[3][i] -
               (REAL)3.0   * zjp1[4][i]);
    zj[4][i] = ONE / (REAL)24.0 * (y[i] -
               (REAL)4.0   * zjp1[1][i] +
               (REAL)6.0   * zjp1[2][i] -
               (REAL)4.0   * zjp1[3][i] +
                             zjp1[4][i]);
  }
  zeigrvek(zj[0], n); zeigrvek(zj[1], n); zeigrvek(zj[2], n);
  zeigrvek(zj[3], n); zeigrvek(zj[4], n);


  /* ------------- adjust step size ------------ */

  for ( ; ; )
  {

    /* - use Newton method for an implicit approximation - */

    for (i = 0; i < n; i++)
      ykp1[i] = zj[0][i] + zj[1][i] + zj[2][i] + zj[3][i] + zj[4][i];

    ykp1[n] = y[n];
    
    dgl(*x + *h, ykp1, f);
    for (k = 0; k < n; k++)
    {
      copy_vector(hilf, ykp1, n);
      hilf[k] -= hs;
      dgl(*x + *h, hilf, fs[k]);
      for (i = 0; i < n; i++)
        fs[k][i] = -*h * (REAL)0.48 * (f[i] - fs[k][i]) / hs;
      fs[k][k] += ONE;
    }
    zeigrmat(fs, n);
    *aufrufe += n + 1;
    for (i = 0; i < n; i++)
    {
      con[i] = ykp1[i] - (REAL)0.48 * (zj[1][i] + TWO * zj[2][i] +
               THREE * zj[3][i] + FOUR * zj[4][i]);
      for (k = 0; k < n; k++)
        fsg[k][i] = fs[i][k];
    }
    zeigrmat(fsg, n);
    fehler = gauss(1, n, fsg, fsg, perm, NULL, NULL, &signdet);
    zeigrmat(fsg, n);
    if (fehler)
    {
      fehler = 6;
      break;
    }
    for (iter = 1; iter <= 3; iter++)
    {
      for (i = 0; i < n; i++)
      {
        hilf[i] = - ykp1[i];
        for (k = 0; k < n; k++)
          hilf[i] += fs[k][i] * ykp1[k];
        hilf[i] = *h * (REAL)0.48 * f[i] + hilf[i] + con[i];
      }
      zeigrvek(hilf, n); zeigrvek(ykp1, n); zeigrmat(fsg, n);
      fehler = gauss(2, n, fsg, fsg, perm, hilf, ykp1, &signdet);
      zeigrmat(fsg, n); zeigrvek(hilf, n); zeigrvek(ykp1, n);
      if (fehler)
      {
        vmfree(vmblock);
        return 8;
      }
      dgl(*x + *h, ykp1, f);
    }
    *aufrufe += 3;

    /* ---- Compute corresponding Gear-Nordsieck approximation ---- */

    for (i = 0; i < n; i++)
      hilf[i] = *h * f[i] - zj[1][i] - TWO * zj[2][i] -
                THREE * zj[3][i] - FOUR * zj[4][i];
    zeigrvek(hilf, n); zeigrvek(zj[4], n);
    for (i = 0; i < n; i++)
    {
      zjp1[0][i] = ykp1[i];
      zjp1[1][i] = *h * f[i];
      zjp1[2][i] = zj[2][i] + THREE * zj[3][i] + SIX * zj[4][i] +
                   (REAL)0.7 * hilf[i];
      zjp1[3][i] = zj[3][i] + FOUR * zj[4][i] + (REAL)0.2 * hilf[i];
      zjp1[4][i] = zj[4][i] + (REAL)0.02 * hilf[i];
    }
    zeigrvek(zjp1[0], n); zeigrvek(zjp1[1], n); zeigrvek(zjp1[2], n);
    zeigrvek(zjp1[3], n); zeigrvek(zjp1[4], n);

    /* -- decide whether to accept last step --*/

    copy_vector(hilf, zjp1[4], n);
    copy_vector(con, zj[4], n);
    zeigrvek(hilf, n); zeigrvek(con, n);
    diff = dist_max(hilf, con, n);
    zeigrvar(diff);
    ymax = norm_max(ykp1, n);
    eps  = (epsabs + epsrel * ymax) / SIX;
    q    = SQRT(SQRT(eps / diff)) / (REAL)1.2;
    if (diff < eps)
    {

      /*  accept last step; prepare for next one --*/

      *x += *h;
      copy_vector(y, ykp1, n);

      /*  stop integration, if interval end xend has been reached  */
      /*  or if there have been too many function calls            */

      nochmal = FALSE;
      do
      {
        if (amEnde)
        {
          *h = dummy;
          vmfree(vmblock);
          return fehler;
        }
        else if (*aufrufe > fmax)
        {
          fehler = 5;
          vmfree(vmblock);
          return fehler;
        }

        /* -- adjust step size for next step --*/

        halt = *h;
        *h   = min(q, TWO) * *h;
        if (*x + *h >= xe)
        {
          dummy  = *h;
          *h     = xend - *x;
          amEnde = 1;

          /* ------- close enough to xend  => stop integration ------ */

          if (*h < eps1 * FABS(xend))
            nochmal = TRUE;
        }
      } while (nochmal);

      /* ------ compute Gear-Nordsieck approximation ---*/
      /* ------ for the next step                    ---*/

      quot1 = *h / halt;
      quot2 = sqr(quot1);
      quot3 = quot2 * quot1;
      quot4 = quot3 * quot1;
      for (i = 0; i < n; i++)
      {
        zj[0][i] = zjp1[0][i];
        zj[1][i] = quot1 * zjp1[1][i];
        zj[2][i] = quot2 * zjp1[2][i];
        zj[3][i] = quot3 * zjp1[3][i];
        zj[4][i] = quot4 * zjp1[4][i];
      }
    }
    else
    {

      /* ------ repeat last step with smaller step size;  ---*/
      /* ------ adjust Gear-Nordsieck approximation   ------ */

      halt  = *h;
      *h    = max(HALF, q) * *h;
      quot1 = *h / halt;
      quot2 = sqr(quot1);
      quot3 = quot2 * quot1;
      quot4 = quot3 * quot1;
      for (i = 0; i < n; i++)
      {
        zj[1][i] = quot1 * zj[1][i];
        zj[2][i] = quot2 * zj[2][i];
        zj[3][i] = quot3 * zj[3][i];
        zj[4][i] = quot4 * zj[4][i];
      }
      amEnde = 0;
    }
  }


  vmfree(vmblock);
  return 0;
}

/* ------------------------------------------------------------------ */

char *gear_fehlertext          /* find error message and error class  */
/*.IX{gear\unt fehlertext}*/
    (
     int       fehlercode,     /* error number from gear()       .....*/
     fehler_t  *fehlerart      /* severity of error from gear()  .....*/
    )                          /* error text .........................*/

/***********************************************************************
* Use error codes from gear() to find the appropriate error message    *
* and classify the type of error as: warning, fatal, no or unknown     *
* error.                                                               *
*                                                                      *
* global variable:                                                     *
* ================                                                     *
* fehler_t                                                             *
***********************************************************************/

{
  static char *fehlertext[] =
    {
      "gear(): no error",
      "gear(): at least one error bound too small",
      "gear(): xend <= x0",
      "gear(): h <= 0",
      "gear(): n <= 0",
      "gear(): too many calls of dgl()",
      "gear(): Jacobi matrix singular",
      "gear(): out of memory",
      "gear(): impossible: 2nd call of gauss() with error"
    };

  switch (fehlercode)
  {
    case 0:  *fehlerart = KEIN_FEHLER;
             break;
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:  *fehlerart = FATAL;
             break;
    default: *fehlerart = UNBEKANNT;
             break;
  }

  return fehlertext[fehlercode];
}

/* -------------------------- END  gear.cpp ------------------------- */
