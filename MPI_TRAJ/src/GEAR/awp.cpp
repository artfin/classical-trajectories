/* ------------------------- MODULE awp.cpp ------------------------- */

/***********************************************************************
*                                                                      *
* Solve an ordinary system of first order differential equations using *
* -------------------------------------------------------------------- *
* automatic step size control                                          *
* ----------------------------                                         *
*                                                                      *
* Programming language: ANSI C                                         *
* Author:               Klaus Niederdrenk (FORTRAN)                    *
* Adaptation:           Juergen Dietel, Computer Center, RWTH Aachen   *
* Source:               existing C, Pascal, QuickBASIC and FORTRAN     *
*                       codes                                          *
* Date:                 6.2.1992, 10.2.1995                            *
*                                                                      *
*                       C++ Release By J-P Moreau, Paris.              *
*                               (www.jpmoreau.fr)                      *
* -------------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges        *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].           *
***********************************************************************/
#include "basis.h"      /*  for  MACH_EPS, FABS, max, SQRT, REAL,     */
                        /*       dglsysfnk, POW, min, FALSE, TRUE,    */
                        /*       boolean, norm_max, ZERO, FOUR, SIX,  */
                        /*       THREE, ONE, TWO, TEN, FIVE, NULL     */
#include "vmblock.h"    /*  for  vmalloc, vmcomplete, vmfree, vminit, */
                        /*       VEKTOR                               */
#include "awp.h"        /*  for  awp, fehler_t, awp_fehlertext        */


/* print a real n vector and its name                                 */
#define zeig(v, n)                 \
  {                                \
    int i;                         \
    printf("%-8s", #v": ");        \
    for (i = 0; i < n; i++)        \
      printf("%16"LZP"g", v[i]);   \
    printf("\n");                  \
  }


/* ------------------------------------------------------------------ */

static REAL dist_max   /* Maximum norm of a difference vector ........*/
        (
         REAL      vektor1[],
         REAL      vektor2[],
         int       n
        )

/***********************************************************************
* Compute the maximum norm of the difference of two [0..n-1] vectors   *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* REAL, FABS, ZERO                                                     *
***********************************************************************/
{
  REAL abstand,        /* reference value for computation of distance */
       hilf;           /* distance of two vector elements             */

  for (n--, abstand = ZERO; n >= 0; n--)
    if ((hilf= FABS(vektor1[n] - vektor2[n])) > abstand)
      abstand = hilf;

  return abstand;
}

/* ------------------------------------------------------------------ */

/***********************************************************************
* Global variables for whole module                                    *
*                                                                      *
* Global variable:                                                     *
* ===============                                                      *
* REAL                                                                 *
***********************************************************************/

static REAL *yhilf,   /* [0..n-1] aux vectors for the embedding       */
            *k1,      /* formulas in ruku23() and engl45().           */
            *k2,      /* dynamically allocated                        */
            *k3,
            *k4,
            *k5,
            *k6;
static REAL *k7,      /* more [0..n-1] aux vectors for embedding      */
            *g6,      /* formula  in prdo45() (dynamic allocation)    */
            *g7;
static int  steif1;   /* Flag, that is set in prdo45() if its         */
                      /* stiffness test (dominant eigenvalue)         */
                      /* indicates so. Otherwise no changes.          */
static int steifanz;  /* counter for number of successive successes   */
                      /* of stiffness test of Shampine and Hiebert in */
                      /* prdo45().                                    */
static int  steif2;   /* Flag, set in prdo45(), when the stiffness    */
                      /* test of  Shampine and Hiebert wa successful  */
                      /* three times in a row; otherwise no changes   */



/* ------------------------------------------------------------------ */

static void ruku23     /* Runge-Kutta embedding f. of 2nd, 3rd degree */
        (
         REAL      x,
         REAL      y[],
         int       n,
         dglsysfnk dgl,
         REAL      h,
         REAL      y2[],
         REAL      y3[]
        )
/***********************************************************************
* Compute 2nd and 3rd order approximates y2, y3 at x + h starting with *
* a solution y at x by using Runge-Kutta embedding formulas on the     *
* first order system of n differential equations   y' = f(x,y) , as    *
* supplied by  dgl().                                                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x    x-value of left end point                                       *
* y    y-value at x                                                    *
* n    number of differential equations                                *
* dgl  pointer to a function that evaluates the right hand side of the *
*      system  y' = f(x,y)                                             *
* h    step size                                                       *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y2   2nd order approximation for y at x + h                          *
* y3   3rd order approximation for y at x + h                          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* yhilf, k1, k2, k3, REAL, dglsysfnk, FOUR, SIX, HALF                  *
***********************************************************************/
{
  int i;                                             /* loop variable */


  (*dgl)(x, y, k1);
  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + h * k1[i];
  (*dgl)(x + h, yhilf, k2);

  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + (REAL)0.25 * h * (k1[i] + k2[i]);
  (*dgl)(x + HALF * h, yhilf, k3);

  for (i = 0; i < n; i++)
    y2[i] = y[i] + HALF * h * (k1[i] + k2[i]),
    y3[i] = y[i] + h / SIX * (k1[i] + k2[i] + FOUR * k3[i]);
}

/* ------------------------------------------------------------------ */

static void engl45     /* Einbettungsforml von England 4. und 5. Ord. */
    (
     REAL      x,                /* Anfangspunkt der Integration .....*/
     REAL      y[],              /* DGLS-Loesung bei x ...............*/
     int       n,                /* Anzahl der DGLen .................*/
     dglsysfnk dgl,              /* rechte Seite des DGLSs ...........*/
     REAL      h,                /* Schrittweite .....................*/
     REAL      y4[],             /* DGLS-Loesung 4. Ordnung bei x+h ..*/
     REAL      y5[]              /* DGLS-Loesung 5. Ordnung bei x+h ..*/
    )

/***********************************************************************
* Compute 4th and 5th order approximates y4, y5 at x + h starting with *
* a solution y at x by using the England embedding formulas on the     *
* first order system of n differential equations   y' = f(x,y) , as    *
* supplied by  dgl().                                                  *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x    initial x-value                                                 *
* y    y-value at x                                                    *
* n    number of differential equations                                *
* dgl  pointer to a function that evaluates the right hand side of the *
*      system  y' = f(x,y)                                             *
* h    step size                                                       *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y4   4th order approximation for y at x + h                          *
* y5   5th order approximation for y at x + h                          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* yhilf, k1, k2, k3, k4, k5, k6, REAL, dglsysfnk, FOUR, SIX, THREE,    *
* TWO, TEN, FIVE                                                       *
***********************************************************************/
{
  int i;                                             /* loop variable */


  (*dgl)(x, y, k1);
  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + HALF * h * k1[i];
  (*dgl)(x + HALF * h, yhilf, k2);

  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + (REAL)0.25 * h * (k1[i] + k2[i]);
  (*dgl)(x + HALF * h, yhilf, k3);

  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + h * (-k2[i] + TWO * k3[i]);
  (*dgl)(x + h, yhilf, k4);

  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + h / (REAL)27.0 * ((REAL)7.0 * k1[i] +
               TEN * k2[i] + k4[i]);
  (*dgl)(x + TWO / THREE * h, yhilf, k5);

  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + h / (REAL)625.0 *
               ((REAL)28.0 * k1[i] - (REAL)125.0 * k2[i] +
                (REAL)546.0 * k3[i] + (REAL)54.0 * k4[i] -
                (REAL)378.0 * k5[i]);
  (*dgl)(x + h / FIVE, yhilf, k6);

  for (i = 0; i < n; i++)
    y4[i] = y[i] + h / SIX * (k1[i] + FOUR * k3[i] + k4[i]),
    y5[i] = y[i] + h / (REAL)336.0 *
            ((REAL)14.0 * k1[i] + (REAL)35.0 * k4[i] +
             (REAL)162.0 * k5[i] + (REAL)125.0 * k6[i]);
}


static void prdo45 /* embedding formulas of Prince-Dormand of 4./5.   */
                   /*  order */
    (
     REAL      x,                /* starting point of integration ....*/
     REAL      y[],              /* initial value at x ...............*/
     int       n,                /* number of DEs    .................*/
     dglsysfnk dgl,              /* right hand side of DE system .....*/
     REAL      h,                /* step size    .....................*/
     REAL      y4[],             /* solution of 4th order at x+h    ..*/
     REAL      y5[]              /* solution of 5th order at x+h    ..*/
    )
/***********************************************************************
* Compute 4th and 5th order approximates y4, y5 at x + h starting with *
* a solution y at x by using the Prince-Dormand embedding formulas on  *
* the first order system of n differential equations y' = f(x,y) , as  *
* supplied by  dgl().                                                  *
* Simultaneously we perform two tests for stiffness whose results are  *
* stored in steif1 and steif2.                                         *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x    initial x-value                                                 *
* y    y-value at x                                                    *
* n    number of differential equations                                *
* dgl  pointer to a function that evaluates the right hand side of the *
*      system  y' = f(x,y)                                             *
* h    step size                                                       *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* y4   4th order approximation for y at x + h                          *
* y5   5th order approximation for y at x + h                          *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* yhilf, k1, k2, k3, k4, k5, k6, k7, g6, g7, REAL, dglsysfnk, THREE,   *
* EIGHT, NINE, steif1, steif2, steifanz                                *
***********************************************************************/
{
  int i;        /* loop variable                                      */
  int steifa;   /* Flag which is set if the second test for stiffness */
                /* (Shampine und Hiebert) is positive; otherwise the  */
                /* flag is erased.                                    */


  (*dgl)(x, y, k1);                           /* coefficients       */
  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + (REAL)0.2 * h *
               k1[i];
  (*dgl)(x + (REAL)0.2 * h, yhilf, k2);

  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + (REAL)0.075 * h *
              (          k1[i]
               + THREE * k2[i]);
  (*dgl)(x + (REAL)0.3 * h, yhilf, k3);

  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + h / (REAL)45.0 *
               (  (REAL) 44.0 * k1[i]
                - (REAL)168.0 * k2[i]
                + (REAL)160.0 * k3[i]);
  (*dgl)(x + (REAL)0.8 * h, yhilf, k4);

  for (i = 0; i < n; i++)
    yhilf[i] = y[i] + h / (REAL)6561.0 *
               (  (REAL)19372.0 * k1[i]
                - (REAL)76080.0 * k2[i]
                + (REAL)64448.0 * k3[i]
                - (REAL) 1908.0 * k4[i]);
  (*dgl)(x + EIGHT / NINE * h, yhilf, k5);

  for (i = 0; i < n; i++)
    g6[i] = y[i] + h / (REAL)167904.0 *
            (  (REAL) 477901.0 * k1[i]
             - (REAL)1806240.0 * k2[i]
             + (REAL)1495424.0 * k3[i]
             + (REAL)  46746.0 * k4[i]
             - (REAL)  45927.0 * k5[i]);
  (*dgl)(x + h, g6, k6);

  for (i = 0; i < n; i++)
    g7[i] = y[i] + h / (REAL)142464.0 *
            (  (REAL)12985.0 * k1[i]
             + (REAL)64000.0 * k3[i]
             + (REAL)92750.0 * k4[i]
             - (REAL)45927.0 * k5[i]
             + (REAL)18656.0 * k6[i]);
  (*dgl)(x + h, g7, k7);

  for (i = 0; i < n; i++)
    y5[i] = g7[i],
    y4[i] = y[i] + h / (REAL)21369600.0 *
            (  (REAL) 1921409.0 * k1[i]
             + (REAL) 9690880.0 * k3[i]
             + (REAL)13122270.0 * k4[i]
             - (REAL) 5802111.0 * k5[i]
             + (REAL) 1902912.0 * k6[i]
             + (REAL)  534240.0 * k7[i]);

  /* Test for stiffness via dominant eigenvalue        */

  if (dist_max(k7, k6, n) > (REAL)3.3 * dist_max(g7, g6, n))
    steif1 = 1;

  /* one step in steffness test of Shampine & Hiebert        */

  for (i = 0; i < n; i++)
    g6[i] = h * ((REAL)2.2   * k2[i] +
                 (REAL)0.13  * k4[i] +
                 (REAL)0.144 * k5[i]),
    g7[i] = h * ((REAL)2.134 * k1[i] +
                 (REAL)0.24  * k3[i] +
                 (REAL)0.1   * k6[i]);
  if (dist_max(g6, g7, n) < dist_max(y4, y5, n))
    steifa = TRUE;
  else
    steifa = FALSE;

  if (steifa)
  {
    steifanz++;
    if (steifanz >= 3)
      steif2 = TRUE;
  }
  else
    steifanz = 0;
}

/* ------------------------------------------------------------------ */
int awp                /* 1st order DESs with autom. step size control*/
        (
         REAL      *x,         /* initial/final x value ..............*/
         REAL      xend,       /* desired end point ..................*/
         int       n,          /* number of DEs ......................*/
         dglsysfnk dgl,        /* right hand side of DEs .............*/
         REAL      y[],        /* initial/final y value ..............*/
         REAL      epsabs,     /* absolute error bound ...............*/
         REAL      epsrel,     /* relative error bound ...............*/
         REAL      *h,         /* initial/final step size ............*/
         int       methode,    /* desired method (3, 6, 7) ...........*/
         long      fmax,       /* maximal # of calls of  dgl() .......*/
         long      *aufrufe    /* actual # of calls of  dgl() ........*/
        )                      /* error code .........................*/

/***********************************************************************
* Compute the solution y of a system of first order ordinary           *
* differential equations       y' = f(x,y)   at xend from the given    *
* initial data (x0, y0).                                               *
* We use automatic step size control internally so that the error of   *
* y (absolutely or relatively) lies within the given error bounds      *
* epsabs and epsrel.                                                   *
*                                                                      *
* Input parameters:                                                    *
* =================                                                    *
* x        initial values                                              *
* y                                                                    *
* n        number of differential equations                            *
* dgl      pointer to a function that evaluates the right hand side    *
*          of the system  y' = f(x,y)                                  *
* xend     end of integration; xend > x0                               *
* h        initial step size                                           *
* epsabs   absolute error bound; >= 0; if = 0 we only check the        *
*          relative error.                                             *
* epsrel   relative error bound; >= 0; if = 0 we check only the        *
*          absolute eror.                                              *
* fmax     max number of evaluations of right hand side in dgl()       *
* methode  chooses the method                                          *
*          = 3: Runge-Kutta method of 2nd/3rd order                    *
*          = 6: England formula of 4th/5th order                       *
*          = 7: Formula of Prince-Dormand of 4th/5th order             *
*                                                                      *
* Output parameters:                                                   *
* ==================                                                   *
* x        final x-value for iteration. If fehler = 0  we usually have *
*            x = xend.                                                 *
* y        final y-value for the solution at x                         *
* h        final step size used; leave for subsequent calls            *
* aufrufe  actual number of calls of dgl()                             *
*                                                                      *
* Return value:                                                        *
* =============                                                        *
* = 0: all ok                                                          *
* = 1: both error bounds chosen too small for the given mach. constant *
* = 2: xend <= x0                                                      *
* = 3: h <= 0                                                          *
* = 4: n <= 0                                                          *
* = 5: more right hand side calls than allowed: aufrufe > fmax,        *
*      x and h contain the current values when stop occured.           *
* = 6: improper input for embedding formula                            *
* = 7: lack of available memory                                        *
* = 8: Computations completed, but the Prince Dormand formula stiff-   *
*      ness test indicates possible stiffness.                         *
* = 9: Computations completed, but both Prince Dormand formula stiff-  *
*      ness tests indicate possible stiffness. Use method for stiff    *
*      systems instead !                                               *
* =10: aufrufe > fmax, see error code 5; AND the Prince Dormand formula*
*      indicates stiffness; retry using a stiff DE solver !            *
*                                                                      *
* global names used:                                                   *
* ==================                                                   *
* ruku23, engl45, dist_max, yhilf, k1, k2, k3, k4, k5, k6, REAL,       *
* MACH_EPS, boolean, FALSE, TRUE, min, max, dglsysfnk, norm_max, FABS, *
* SQRT, vminit, vmalloc, vmcomplete, vmfree, VEKTOR, POW, ZERO, ONE,   *
* TWO, NULL, prdo45, steif1, steif2, steifanz                          *
***********************************************************************/
{
#define MACH_2   ((REAL)100.0 * MACH_EPS)  /* machine constant related*/
                                           /* value used for break-off*/
                                           /* criteria as zero        */

  REAL    xend_h,    /* |xend| - MACH_2, carrying same sign as xend   */
          ymax,      /* Maximum norm of newest approximation of max   */
                     /* order                                         */
          hhilf = ZERO,  /* aux storage for the latest value of h     */
                     /* produced by step size control. It is saved    */
                     /* here in order to avoid to return a `h' that   */
                     /* resulted from an arbitrary reduction at the   */
                     /* end of the interval.                          */
          diff,      /* distance of the two approximations from the   */
                     /* embedding formula                             */
          s,         /* indicates acceptance level for results from   */
                     /* embeding formula                              */
          *y_bad,    /* approximate solution of low order             */
          *y_good,   /* ditto of high order                           */
          mach_1;    /* machine constant dependent variable which     */
                     /* avoids using too little steps near xend       */
  int     i,         /* Loop variable                                 */
          fehler;    /* error code                                    */
  boolean amEnde,    /* flag that shows if the end of the interval    */
                     /* can be reached with the actual step size      */
          fertig;    /* Flag indicating end of iterations             */
  void    *vmblock,  /* List of dynamic allocations                   */
          (*einbett)(REAL      x,      /* pointer to the function     */
                     REAL      *y,     /* that computes the two       */
                     int       n,      /* approximations via a pair   */
                     dglsysfnk dgl,    /* of embedding formulas       */
                     REAL      h,
                     REAL      *ybad,
                     REAL      *ygut) = NULL;


  fehler   = 0;                         /* initialize some variables  */
  mach_1   = POW(MACH_EPS, (REAL)0.75);
  amEnde   = FALSE;
  fertig   = FALSE;
  steif1   = FALSE;
  steif2   = FALSE;
  steifanz = 0;
  *aufrufe = 0l;
  ymax     = norm_max(y, n);
  xend_h   = (xend >= ZERO) ? (xend * (ONE - MACH_2)) :
                              (xend * (ONE + MACH_2));

  switch (methode)            /* point function to correct embedding  */
  {                           /* formula                              */
    case 3:
      einbett = ruku23;
      break;
    case 6:
      einbett = engl45;
      break;
    case 7:
      einbett = prdo45;
  }


  /* ------------------- check input               ------------------ */

  if (epsabs <= MACH_2 * ymax && epsrel <= MACH_2)
    return 1;
  if (xend_h < *x)
    return 2;
  if (*h < MACH_2 * FABS(*x))
    return 3;
  if (n <= 0)
    return 4;
  if (methode != 3 && methode != 6 && methode != 7)
    return 6;


  /* - allocate aux vectors : six [0..n-1] vectors for method = 3,  - */
  /* - nine for  methode = 6 or 7                                   - */

  vmblock = vminit();                 /* initialize storage           */
  y_bad  = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  y_good = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  yhilf  = (REAL *)vmalloc(vmblock, VEKTOR, n+1, 0);
  yhilf[n] = y[n];
  k1     = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  k2     = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  k3     = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (methode == 6 || methode == 7)
    k4   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0),
    k5   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0),
    k6   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (methode == 7)
    k7   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0),
    g6   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0),
    g7   = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))                      /* lack of memory ? */
  {
    vmfree(vmblock);                /* free storage and report error  */
    return 7;
  }

  /*********************************************************************
  *                                                                    *
  *                 I t e r a t i o n                                  *
  *                                                                    *
  *********************************************************************/

  if (*x + *h > xend_h)            /* almost at end point ?          */
    hhilf  = *h,                    /* A shortened step might be      */
    *h     = xend - *x,            /* enough.                        */
    amEnde = TRUE;

  do                       /* solve DE system by integrating from     */
  {                        /* x0 to xend by suitable steps            */

    (*einbett)(*x, y, n, dgl, *h, y_bad, y_good);     /* integrate   */
#ifdef DEBUG
    zeig(y_bad, n);
    zeig(y_good, n);
#endif
    *aufrufe += methode;

    if ((diff = dist_max(y_bad, y_good, n)) < MACH_2)  /* compute s   */
      s = TWO;
    else
    {
      ymax = norm_max(y_good, n);
      s    = SQRT(*h * (epsabs + epsrel * ymax) / diff);
      if (methode != 3)
        s = SQRT(s);
    }

    if (s > ONE)               /* integration acceptable? */
    {
      for (i = 0; i < n; i++)  /* accept highest order solution       */
        y[i] = y_good[i];      /* move x                             */
      *x += *h;

      if (amEnde)                /* at end of interval?               */
      {
        fertig = TRUE;           /* stop iteration                    */
        if (methode == 7)
        {
          if (steif1 || steif2)
            fehler = 8;
          if (steif1 && steif2)
            fehler = 9;
        }
      }
      else if (*aufrufe > fmax)  /* too many calls of   dgl()?        */
      {
        hhilf  = *h,             /* save actual step size             */
        fehler = 5,              /* report error and stop             */
        fertig = TRUE;
        if (methode == 7 &&
            (steif1 || steif2))
          fehler = 10;
      }

      else                          /* Integration was successful     */
      {                             /* not at the interval end?       */
        *h *= min(TWO,              /* increase step size for next    */
                  (REAL)0.98 * s);  /* step properly, at most by      */
                                    /* factor two. Value `0.98*s' is  */
                                    /* used in order to avoid that    */
                                    /* the theoretical value s is     */
                                    /* exceeded by accidental         */
                                    /* rounding errors.               */
        if (*x + *h > xend_h)      /* nearly reached xend?           */
        {
          hhilf  = *h;              /* => One further step with       */
          *h     = xend - *x;      /*    reduced step size might be  */
          amEnde = TRUE;            /*    enough.                     */

          if (*h < mach_1 * FABS(xend))    /* very close to xend ?    */
            fertig = TRUE;                 /* finish iteration        */
        }
      }
    }

    else                          /* step unsuccessful?               */
    {                             /* before repeating this step:      */
      *h *= max(HALF,             /* reduce step size properly, at    */
                (REAL)0.98 * s);  /* most by factor 1/2 (for factor   */
      amEnde = FALSE;             /* 0.98: see above)                 */
    }

  }
  while (! fertig);


  vmfree(vmblock);
  *h = hhilf;      /* return the latest step size computed by step    */
                   /* size control and                                */
  return fehler;   /* and error code to the caller                    */
}


/* ------------------------------------------------------------------ */

char *awp_fehlertext   /* determine error code and class              */
    (
     int       fehlercode,     /* Number of error from awp()   .......*/
     fehler_t  *fehlerart      /* type of error from awp()      ......*/
    )                          /* error text .........................*/
/***********************************************************************
* Print out the appropriate error text and its severity fro mawp().    *
* Sample output: All ok, Caution, fatal, unknown error.                *
*                                                                      *
* Global name:                                                         *
* ============                                                         *
* fehler_t                                                             *
***********************************************************************/
{
  static char *fehlertext[] =
    {
      "awp(): no error",
      "awp(): error bounds too small",
      "awp(): xend <= x0",
      "awp(): h <= 0",
      "awp(): n <= 0",
      "awp(): too many calls of dgl()",
      "awp(): improper embedding formula",
      "awp(): out of memory",
      "awp(): caution: stiffness suspected (one test)",
      "awp(): caution: stiffness suspected (two tests)",
      "awp(): too many calls of dgl() and stiffness suspected"
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
    case 10: *fehlerart = FATAL;
             break;
    case 8:
    case 9:  *fehlerart = WARNUNG;
             break;
    default: *fehlerart = UNBEKANNT;
             break;
  }

  return fehlertext[fehlercode];
}

/* -------------------------- END  awp.cpp -------------------------- */
