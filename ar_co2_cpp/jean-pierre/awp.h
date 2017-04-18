/*.Adaptive Methods for Initial Value Problems [See BIBLI 11]
   Automatic Step Size Control,
   Adaptive Methods for Initial Value Problems */

/* ------------------------ DECLARATIONS awp.h ---------------------- */

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
        );                     /* error code .........................*/

char *awp_fehlertext          /* determine error type and its gravity */
    (
     int       fehlercode,     /* error number from awp() ............*/
     fehler_t  *fehlerart      /* gravity of error from  awp() .......*/
    );                         /* error message ......................*/


/* ---------------------------- END awp.h --------------------------- */
