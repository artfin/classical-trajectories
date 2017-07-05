

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




#include "matrix.h"





void 
syst (REAL t ,REAL *y, REAL *f)
{
  (void)(t); // avoid unused parameter warning 
  //double * par = (double *)params;

  

  double * out = new double[6];
  rhs(out, y[0],y[1],y[2],y[3],y[4],y[5], y[6]);// R  Theta pR pT phi theta
 //hamiltonian(out, y[0],y[2],y[1],y[3],y[4],y[5],10.0);// R pR Theta pT phi theta := y[0]..y[5]

  f[0] = out[0];//dR/dt  
 // f[1]=f1(y,par)-(dpotdR);//d(pR)/dt 
  f[1]=out[1]; //d(Theta)/dt
 
  f[2] = out[2]; //d(pR)/dt
  f[3] = out[3];//d(pT)/dt
 
  f[4] = out[4];//d(phi)/dt
  f[5] = out[5];//d(theta)/dt



  delete [] out;
  //delete [] par;

  
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

  FILE * trajectory_file = fopen("traj_cpp_J.txt","w");

/* -------------------- read input  -------------------- */

  N = 6;
 
  vmblock = vminit();                 /* initialize storage */
  y0  = (REAL *)vmalloc(vmblock, VEKTOR, N+1, 0);
 
  if (! vmcomplete(vmblock))  {       /* out of memory? */
    printf("mgear: out of memory.\n");
    return 0;
  }

  double step = 50000;
  double end = 5000000;

  int num_points = (end/step);




  epsabs = 1E-13;
  epsrel = 1E-13;
  

  t0 = 0.0;

   
  y0[0] =15.0;
  y0[1] = 1.0;
  y0[2] = -4.0;
  y0[3] = 1.0;
  y0[4] = 1.0;
  y0[5] = 1.0;

  y0[6] = 1.0;



  h =0.1;
  xend = step;
  fmax = 1000000; 



  for (int i = 0; i < num_points; ++i)
  {
     fehler = gear4(&t0, xend, N, syst, y0, epsabs,
                 epsrel, &h, fmax, &aufrufe);
     if (fehler != 0) {
     printf(" Gear4: error n° %d\n", 10 + fehler);
     return 0;
     }

     
   // if ( i == 0)
     // fprintf(trajectory_file,"%f %.12f %.12f %.12f %.12f %.12f\n",
       //     t0, y0[0],y0[1],y0[2],y0[3],y0[4],y0[5]);
     fprintf(trajectory_file,"%f %.12f %d\n", t0, y0[0], aufrufe);
      xend = step*(i+2);
      aufrufe = 0;

  

  }
 

  

  for (i = 0; i < 1; i++)
  {
   // if ( i == 0)
      printf("approximate solution y%d(x) = %24.15"LZP"e\n",
            i + 1, y0[i]);

  }

  
  fclose(trajectory_file);
  return 0;
}

/* -------------------------- END mgear.cpp ------------------------- */
