//* ------------------------ MODULE fgauss.cpp ----------------------- */

/****************************************************************
*         Gauss algorithm for solving linear equations          *
* ------------------------------------------------------------- *
* REF.: "Numerical Algorithms with C, By Gisela Engeln-Muellges *
*        and Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].    *
*                                                               *
*                           C++ Release By J-P Moreau, Paris.   *
*                                   (www.jpmoreau.fr)           *
****************************************************************/ 
#include "basis.h"
#include "vmblock.h"

//headers of functions used by gauss()
int gaudec              /* Gauss decomposition .......................*/
           (
            int     n,            /* size of matrix ..................*/
            REAL *  mat[],        /* Input matrix ....................*/
            REAL *  lumat[],      /* matrix decomposition ............*/
            int     perm[],       /* row interchanges ................*/
            int *   signd         /* sign of perm ....................*/
           );

int gausol              /* Gauss solution ............................*/
           (
            int     n,            /* size of matrix ..................*/
            REAL *  lumat[],      /* decomposed matrix (LU) ..........*/
            int     perm[],       /* row permutation vector ..........*/
            REAL    b[],          /* Right hand side .................*/
            REAL    x[]           /* solution ........................*/
           );


int gauss            /* Gauss algorithm for solving linear equations .*/
          (
           int     mod,           /* Modus: 0, 1, 2, 3 ...............*/
           int     n,             /* Dimension of matrix .............*/
           REAL *  mat[],         /* Input matrix ....................*/
           REAL *  lumat[],       /* LU decomposition ................*/
           int     perm[],        /* row remutation vector ...........*/
           REAL    b[],           /* right hand side .................*/
           REAL    x[],           /* solution of the system ..........*/
           int *   signd          /* sign of the permutation .........*/
          )
/*====================================================================*
 *                                                                    *
 *  The funktion gauss solves a linear system :  mat * x = b.         *
 *  Here mat is the nonsingular system matrix, b the right hand side  *
 *  of the system and x the solution vector.                          *
 *                                                                    *
 *  gauss uses the Gauss algorithm and computes a triangular factori- *
 *  zation of mat and skaled column pivot search.  (Crout method with *
 *  row swaps).                                                       *
 *                                                                    *
 *====================================================================*
.BE*)
 *                                                                    *
 *   Application:                                                     *
 *   ============                                                     *
 *      Solve general linear system with a nonsingular coefficient    *
 *      matrix.                                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Control parameter:                                               *
 *   ==================                                               *
 *      mod      int mod;                                             *
 *               calling modus for gauss:                             *
 *       = 0     Find factorization and solve linear system           *
 *       = 1     Find factorization only.                             *
 *       = 2     Solve linear system only; the factorization is       *
 *               already available in lumat. This saves work when     *
 *               solving a linear system repeatedly for several right *
 *               hand sides and the same system matrix such as when   *
 *               inverting the matrix.                                *
 *       = 3     as under 2, additionally we improve the solution     *
 *               via iterative refinement.                            *
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of mat and lumat,                          *
 *               size of the vector b, the right hand side, the       *
 *               solution x and the permutation vector perm.          *
 *      mat      REAL   *mat[];                                       *
 *               Matrix of the linear system. It is stored in vector  *
 *               form.                                                *
 *      lumat    REAL   *lumat[];      ( for mod = 2, 3 )             *
 *               LU factors of mat                                    *
 *               mat in eine untere und obere Dreieckmatrix ent-      *
 *               lumat can be stored in the space of mat.             *
 *      perm     int perm[];           ( for mod = 2, 3 )             *
 *               Permutation vector, of the row interchangfes in lumat*
 *      b        REAL   b[];           ( for mod = 0, 2, 3 )          *
 *               Right hand side of the system.                       *
 *      signd    int *signd;           ( for mod = 2, 3 )             *
 *               sign of the permutation in perm; the determinant of  *
 *               mat can be computed as the product of the diagonal   *
 *               entries of lumat times signd.                        *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      lumat    REAL   *lumat[];      ( for mod = 0, 1 )             *
 *               LU factorization of mat.                             *
 *      perm     int perm[];           ( for mod = 0, 1 )             *
 *               row ermutation vektor                                *
 *      x        REAL   x[];           ( for mod = 0, 2, 3 )          *
 *               solution vector.                                     *
 *      signd    int *signd;           ( for mod = 0, 1 )             *
 *               sign of perm.                                        *
 *                                                                    *
 *   Return value :                                                   *
 *   ==============                                                   *
 *      =-1      Max. number (MAXITER) of iterative refinements       *
 *               reached (MAXITER) while mod = 3                      *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or other invalid input                         *
 *      = 2      lack of memory                                       *
 *      = 3      Matrix singular                                      *
 *      = 4      Matrix numerically singular                          *
 *      = 5      incorrect call                                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions used :                                                 *
 *   ================                                                 *
 *                                                                    *
 *      int gaudec  (): determines  LU decomposition                  *
 *      int gausol  (): solves the linear system                      *
 *      int gausoli (): solves the system with iterative refinement   *
 *                                                                    *
 *====================================================================*/
{
  int  rc;

  if (n < 1) return (1);

  switch (mod)
  {
    case 0: /* Find factorization and solve system ...................*/
            rc = gaudec (n, mat, lumat, perm, signd);
            if (rc == 0)
              return (gausol (n, lumat, perm, b, x));
            else
              return (rc);

    case 1: /* Find factorization only ...............................*/
            return (gaudec (n, mat, lumat, perm, signd));

    case 2: /* Solve only ............................................*/
            return (gausol (n, lumat, perm, b, x));

    case 3: /* solve and then use iterative refinement ...............*/
            //return (gausoli (n, mat, lumat, perm, b, x));
		    printf(" fgauss: gausoli not implemented.\n");
			return 0;
  }

  return (5);                                           /* Wrong call */
}

int gaudec              /* Gauss decomposition .......................*/
           (
            int     n,            /* size of matrix ..................*/
            REAL *  mat[],        /* Input matrix ....................*/
            REAL *  lumat[],      /* matrix decomposition ............*/
            int     perm[],       /* row interchanges ................*/
            int *   signd         /* sign of perm ....................*/
           )
/*====================================================================*
 *                                                                    *
 *  gaudec decomposes a nonsingular n x n matrix into a product of a  *
 *  unit lower and an upper triangular matrix. Both triangular factors*
 *  are stored in lumat (minus the unit diagonal, of course).         *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Eingabeparameter:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of  mat and lumat,                         *
 *               size of  b , x and perm                              *
 *      mat      REAL   *mat[];                                       *
 *               original system matrix in vector form.               *
 *                                                                    *
 *   Output parameters:                                               *
 *   ==================                                               *
 *      lumat    REAL   *lumat[];                                     *
 *               LU factorization                                     *
 *      perm     int perm[];                                          *
 *               row permutation vector for lumat                     *
 *      signd    int *signd;                                          *
 *               sign of perm. The determinant of mat can be computed *
 *               as the product of the diagonal entreis of lumat times*
 *               signd.                                               *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or invalid input                               *
 *      = 2      lack of memory                                       *
 *      = 3      Matrix is singular                                   *
 *      = 4      Matrix numerically singular                          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use :                                               *
 *   ==================                                               *
 *                                                                    *
 *      void *vmalloc():  allocate vectors or matrices                *
 *      void vmfree():    free list of vectors and matrices           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants in use  :  NULL, MACH_EPS                              *
 *   ===================                                              *
 *                                                                    *
 *   Macros: SWAP, ABS                                                *
 *   ======                                                           *
 *====================================================================*/
{
  int  m, j, i, j0;
  REAL piv, tmp, *d = NULL, zmax;
  void *vmblock;

  if (n < 1) return (1);                   /*  Invalid parameters     */

  if (mat == NULL || lumat == NULL) return (1);
  if (perm == NULL) return (1);

  for (i = 0; i < n; i++)
    if (mat[i] == NULL || lumat[i] == NULL) return (1);

                                  /* d = Skaling vector for pivoting  */
  vmblock = vminit();                   /* allocate storage           */
  d = (REAL *)vmalloc(vmblock, VEKTOR, n, 0);
  if (! vmcomplete(vmblock))            /* lack of memory?            */
    return 2;

  if (lumat != mat)                     /* If  lumat and  mat are     */
    CopyMat (n, n, mat, lumat);         /* distinct, copy mat to lumat*/

  for (i = 0; i < n; i++)
  {
    perm[i] = i;                        /* Initialize perm            */
    zmax = ZERO;
    for (j = 0; j < n; j++)             /* find row maxima            */
    {
      tmp = ABS (lumat[i][j]);
      if (tmp > zmax) zmax = tmp;
    }

    if (zmax == ZERO)                   /* mat is singular            */
    {
      vmfree(vmblock);
      return (3);
    }
    d[i] = ONE / zmax;
  }

  *signd = 1;                         /* initialize sign of perm      */

  for (i = 0; i < n; i++)
  {
    piv = ABS (lumat[i][i]) * d[i];
    j0 = i;                           /* Search for pivot element     */
    for (j = i + 1; j < n; j++)
    {
      tmp = ABS (lumat[j][i]) * d[j];
      if (piv < tmp)
      {
        piv = tmp;                    /* Mark pivot element and       */
        j0 = j;                       /* its location                 */
      }
    }

    if (piv < MACH_EPS)               /* If piv is small, mat is      */
    {                                 /* nearly singular              */
      *signd = 0;
      vmfree(vmblock);
      return (4);
    }

    if (j0 != i)
    {
      *signd = - *signd;              /* update signd                 */

      SWAP (int, perm[j0], perm[i]);  /* swap pivotentries            */
      SWAP (REAL, d[j0], d[i])        /* swap skaling vector          */

      SWAP (REAL *, lumat[j0], lumat[i]);    /* swap j0-th and i-th   */
                                             /* row of  lumat         */
    }

    for (j = i + 1; j < n; j++)       /* Gauss elimination            */
    {
      if (lumat[j][i] != ZERO)
      {
        lumat[j][i] /= lumat[i][i];
        tmp = lumat[j][i];
        for (m = i + 1; m < n; m++)
          lumat[j][m] -= tmp * lumat[i][m];
      }
    }
  } /* end i */

  vmfree(vmblock);                   /* free space of scaling vector  */
  return (0);
}

int gausol              /* Gauss solution ............................*/
           (
            int     n,            /* size of matrix ..................*/
            REAL *  lumat[],      /* decomposed matrix (LU) ..........*/
            int     perm[],       /* row permutation vector ..........*/
            REAL    b[],          /* Right hand side .................*/
            REAL    x[]           /* solution ........................*/
           )
/*====================================================================*
 *                                                                    *
 *  gausol  finds the solution x of the linear system  lumat * x = b  *
 *  for the product matrix lumat, that describes an LU decomposition, *
 *  as produced by gaudec.                                            *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of lumat,                                  *
 *      lumat    REAL   *lumat[];                                     *
 *               LU factorization, as produced from gaudec            *
 *      perm     int perm[];                                          *
 *               row permutation vector for lumat                     *
 *      b        REAL   b[];                                          *
 *               Right hand side of the system.                       *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      x        REAL   x[];                                          *
 *               solution vector                                      *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or other invalid input parameter               *
 *      = 3      improper LU decomposition ( zero diagonal entry)     *
 *                                                                    *
 *====================================================================*/
{
  int  j, k;
  REAL sum;

  if (n < 1) return (1);                   /* Invalid input parameter */

  if (lumat == NULL || b == NULL || perm == NULL) return (1);

  for (j = 0; j < n; j++)
    if (lumat[j] == NULL) return (1);

  for (k = 0; k < n; k++)                              /* update b    */
  {
    x[k] = b[perm[k]];
    for (j = 0; j < k; j++)
      x[k] -= lumat[k][j] * x[j];
  }

  for (k = n - 1; k >= 0; k--)                    /* back substitute  */
  {
    sum = ZERO;
    for (j = k + 1; j < n; j++)
      sum += lumat[k][j] * x[j];

    if (lumat[k][k] == ZERO) return (3);
    x[k] = (x[k] - sum) / lumat[k][k];
  }

  return (0);
}

REAL det                /* Determinant  ..............................*/
           (
            int     n,            /* Dimension of the matrix .........*/
            REAL *  mat[]         /* matrix ..........................*/
           )
/*====================================================================*
 *                                                                    *
 *  det computes the determinant of an n x n real matrix mat          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameter:                                                 *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of mat                                     *
 *      mat      REAL   *mat[];                                       *
 *               n x n matrix                                         *
 *                                                                    *
 *   Return value:                                                    *
 *   =============                                                    *
 *      REAL     Determinant of mat.                                  *
 *               If the return value = 0, then the matrix is singular *
 *               or the storage is insufficient                       *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Functions in use  :                                              *
 *   ===================                                              *
 *                                                                    *
 *      int gaudec ():    LU decomposition of mat                     *
 *      void *vmalloc():  allocate vectors or matrices                *
 *      void vmfree():    free list of vectors and matrices           *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Constants in use  :  NULL, MAXROOT, EPSQUAD                      *
 *   ===================                                              *
 *                                                                    *
 *====================================================================*/
{
  int   i, rc, signd, *perm;
  REAL  **lu, tmpdet;
  void *vmblock;

  if (n < 1) return (ZERO);
                                              /* create buffer for    */
  vmblock = vminit();                         /* LU decomposition     */
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));

  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return (ZERO);
  }

  rc = gaudec (n, mat, lu, perm, &signd);     /* decompose            */

  if (rc != 0 || signd == 0)
  {
    vmfree(vmblock);
    return (ZERO);
  }

  tmpdet = (REAL) signd;
  for (i = 0; i < n; i++)
  {
    if (ABS(tmpdet) < EPSQUAD)
    {
      vmfree(vmblock);
      return (ZERO);
    }
    else
    if (ABS(tmpdet) > MAXROOT || ABS(lu[i][i]) > MAXROOT)
    {
      vmfree(vmblock);
      return (MAXROOT);
    }
    else
      tmpdet *= lu[i][i];                     /* compute det          */
  }

  vmfree(vmblock);

  return (tmpdet);
}

int mgauss                     /* Gauss for multiple right hand sides */
           (
            int     n,            /* Dimension of system .............*/
            int     k,            /* number of right hand sides ......*/
            REAL *  mat[],        /* original matrix .................*/
            REAL *  rmat[]        /* Right hand sides/solutions ......*/
           )
/*====================================================================*
 *                                                                    *
 *  mgauss  finds the solution matrix x for the linear system         *
 *  mat * x = rmat with an  n x n coefficient matrix mat and a        *
 *  n x k matrix rmat of right hand sides. Here mat must be           *
 *  nonsingular.                                                      *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   Input parameters:                                                *
 *   ================                                                 *
 *      n        int n;  ( n > 0 )                                    *
 *               Dimension of mat.                                    *
 *      k        int k;  ( k > 0 )                                    *
 *               number of right hand sides                           *
 *      mat      REAL   *mat[];                                       *
 *               n x n original system matrix                         *
 *      rmat     REAL   *rmat[];                                      *
 *               matrix of right hand sides                           *
 *                                                                    *
 *   Output parameter:                                                *
 *   ================                                                 *
 *      rmat     REAL   *rmat[];                                      *
 *               solution matrix for the system.                      *
 *               The input right hand sides are lost.                 *
 *                                                                    *
 *   Return value :                                                   *
 *   =============                                                    *
 *      = 0      all ok                                               *
 *      = 1      n < 1 or k < 1 or invalid input parameter            *
 *      = 2      lack of memory                                       *
 *      = 3      mat is numerically singular                          *
 *                                                                    *
 *====================================================================*
 *                                                                    *
 *   functions used :                                                 *
 *   ================                                                 *
 *                                                                    *
 *      int gaudec ():    LU decomposition of mat                     *
 *      void *vmalloc():  allocate vector or matrix                   *
 *      void vmfree():    free list of vectors or matrix              *
 *                                                                    *
 *====================================================================*/
{
  register int i, j;
  int      m, *perm, signd, rc;
  REAL     **lu, *x, sum;
  void *vmblock;

  if (n < 1 || k < 1) return (1);                /* Invalid parameter */

  if (mat == NULL || rmat == NULL) return (1);
  if (mat == rmat) return (1);

  for (j = 0; j < n; j++)
    if (mat[j] == NULL || rmat[j] == NULL) return (1);
                                              /* allocate storage     */
  vmblock = vminit();                         /* for LU factorization */
  lu   = (REAL **)vmalloc(vmblock, MATRIX,  n, n);
  perm = (int *)  vmalloc(vmblock, VVEKTOR, n, sizeof(*perm));
  x    = (REAL *) vmalloc(vmblock, VEKTOR,  n, 0);

  if (! vmcomplete(vmblock))
  {
    vmfree(vmblock);
    return (2);
  }

  rc = gaudec (n, mat, lu, perm, &signd);     /* compute factorization*/
                                              /*  in lu               */
  if (rc != 0 || signd == 0)                  /* if not possible      */
  {
    vmfree(vmblock);                          /* release storage      */
    return (3);
  }

  for (m = 0; m < k; m++)          /* Loop over the right hand sides  */
  {
    for (i = 0; i < n; i++)                 /* Updating the b's       */
    {
      x[i] = rmat[perm[i]][m];
      for (j = 0; j < i; j++)
        x[i] -= lu[i][j] * x[j];
    }

    for (i = n - 1; i >= 0; i--)            /* back substitution      */
    {
      sum = ZERO;
      for (j = i + 1; j < n; j++)
        sum += lu[i][j] * x[j];

      if (lu[i][i] == ZERO)               /* invalid LU decomposition */
      {
        rc = 2;
        break;
      }
      x[i] = (x[i] - sum) / lu[i][i];
    }

    for (j = 0; j < n; j++)                 /* Save result            */
      rmat[j][m] = x[j];
  }

  vmfree(vmblock);                          /* free storage           */

  return (rc);
}

/* ------------------------- END fgauss.cpp ------------------------- */