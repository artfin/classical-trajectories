/* ---------------------- DECLARATIONS t_dgls.h --------------------- */

#ifndef T_DGLS_H_INCLUDED
#define T_DGLS_H_INCLUDED

typedef struct
{
  int       n;                             /* number of equations,    */
  dglsysfnk rechte_seite;                  /* right hand side,        */
  char      *(*dgl_text)(void);            /* right hand side as text */
  void      (*exakte_loesung)(REAL x,      /* and analytic solution   */
                              REAL y[]     /* of the DE system        */
                             );
} bsptyp;

bsptyp *dgls_waehlen(int nummer);

#endif

/* -------------------------- END t_dgls.h -------------------------- */
