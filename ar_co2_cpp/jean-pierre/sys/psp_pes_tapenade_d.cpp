/*        Generated by TAPENADE     (INRIA, Ecuador team)
    Tapenade 3.12 (r6213) - 13 Oct 2016 10:54
*/
#include "psp_pes_tapenade_d.h"

/*
  Differentiation of potential in forward (tangent) mode:
   variations   of useful results: *result
   with respect to varying inputs: Theta R
   RW status of diff variables: *result:out Theta:in R:in
   Plus diff mem management of: result:in
*/
void potential_d(double R, double Rd, double Theta, double Thetad, double &
        result, double &resultd) {
    double V0, V2, V4, V6, V8, V10;
    double V0d, V2d, V4d, V6d, V8d, V10d;
    double arg1;
    double arg1d;
    //double R = vars[0];
    //ouble Theta = vars[1];
    arg1d = -(0.716288*Rd) - 0.0869136*(Rd*R+R*Rd);
    arg1 = -0.716288*R - 0.0869136*R*R;
    V0d = 22.4247*arg1d*exp(arg1);
    V0 = 22.4247*exp(arg1);
    if (R >= 6.32925) {
        V0d = V0d - -(2380*Rd*8*pow(R, 8-1))/(pow(R, 8)*pow(R, 8)) + 114.5*Rd*
            6*pow(R, 6-1)/(pow(R, 6)*pow(R, 6));
        V0 -= 114.5/pow(R, 6) + 2380/pow(R, 8);
    } else {
        arg1d = -(0.650479*Rd) - 0.0320299*(Rd*R+R*Rd);
        arg1 = -0.650479*R - 0.0320299*R*R;
        V0d = V0d - 0.599056*arg1d*exp(arg1);
        V0 += -0.599056*exp(-0.650479*R-0.0320299*R*R);
    }
    arg1d = -(0.811806*Rd) - 0.075313*(Rd*R+R*Rd);
    arg1 = -0.811806*R - 0.075313*R*R;
    V2d = 63.5744*arg1d*exp(arg1);
    V2 = 63.5744*exp(arg1);
    if (R < 6.90026) {
        arg1d = -(0.620877*Rd) - 0.0310717*(Rd*R+R*Rd);
        arg1 = -0.620877*R - 0.0310717*R*R;
        V2d = V2d - 0.207391*arg1d*exp(arg1);
        V2 -= 0.207391*exp(-0.620877*R-0.0310717*R*R);
    } else {
        V2d = V2d - -(2080*Rd*8*pow(R, 8-1))/(pow(R, 8)*pow(R, 8)) + 26.6*Rd*6
            *pow(R, 6-1)/(pow(R, 6)*pow(R, 6));
        V2 -= 26.6/pow(R, 6) + 2080/pow(R, 8);
    }
    arg1d = -(1.17577*Rd) - 0.0477138*(Rd*R+R*Rd);
    arg1 = -1.17577*R - 0.0477138*R*R;
    V4d = 99.1128*arg1d*exp(arg1);
    V4 = 99.1128*exp(arg1);
    if (R < 6.96450) {
        arg1d = -(0.73534*Rd) - 0.0296750*(Rd*R+R*Rd);
        arg1 = -0.73534*R - 0.0296750*R*R;
        V4d = V4d - 0.0523497*arg1d*exp(arg1);
        V4 -= 0.0523497*exp(-0.73534*R-0.0296750*R*R);
    } else {
        V4d = V4d - -(410*Rd*8*pow(R, 8-1))/(pow(R, 8)*pow(R, 8));
        V4 -= 410/pow(R, 8);
    }
    arg1d = -(1.05547*Rd) - 0.0182219*(Rd*R+R*Rd);
    arg1 = -1.05547*R - 0.0182219*R*R;
    V6d = -(318.652*1.88135*Rd*exp(-1.88135*R)) - 0.0374994*arg1d*exp(arg1);
    V6 = 318.652*exp(-1.88135*R) - 0.0374994*exp(arg1);
    arg1d = -(1.03942*Rd) - 0.0494781*(Rd*R+R*Rd);
    arg1 = -1.03942*R - 0.0494781*R*R;
    V8d = -(332.826*2.14596*Rd*exp(-2.14596*R)) - 0.0137376*arg1d*exp(arg1);
    V8 = 332.826*exp(-2.14596*R) - 0.0137376*exp(arg1);
    V10d = 0.108283*2.04765*Rd*exp(-2.04765*R) - 435.837*2.44616*Rd*exp(-
        2.44616*R);
    V10 = 435.837*exp(-2.44616*R) - 0.108283*exp(-2.04765*R);
    double cosine;
    double cosined;
    cosined = -(Thetad*sin(Theta));
    cosine = cos(Theta);
    double pes;
    double cosine2 = cosine*cosine;
    double cosine2d = cosined*cosine + cosine*cosined;
    double cosine4 = cosine2*cosine2;
    double cosine4d = cosine2d*cosine2 + cosine2*cosine2d;
    double cosine6 = cosine2*cosine*(cosine2*cosine);
    double cosine6d = (cosine2d*cosine+cosine2*cosined)*cosine2*cosine + 
    cosine2*cosine*(cosine2d*cosine+cosine2*cosined);
    double cosine10 = cosine4*cosine*(cosine4*cosine);
    double cosine10d = (cosine4d*cosine+cosine4*cosined)*cosine4*cosine + 
    cosine4*cosine*(cosine4d*cosine+cosine4*cosined);
    /*pes = gsl_sf_legendre_Pl(0,cosine)*V0 + gsl_sf_legendre_Pl(2,cosine)*V2 + gsl_sf_legendre_Pl(4,cosine)*V4 + 
			  gsl_sf_legendre_Pl(6,cosine)*V6 + gsl_sf_legendre_Pl(8,cosine)*V8 + gsl_sf_legendre_Pl(10,cosine)*V10;  
    */
    resultd = V0d + 3.0*cosine2d*V2/2 + (3.0/2*cosine2-1.0/2)*V2d + (35.0*
        cosine4d/8-15.0*cosine2d/4)*V4 + (35.0/8*cosine4-15.0/4*cosine2+3.0/8)
        *V4d + (231.0*cosine6d/16-315.0*cosine4d/16+105.0*cosine2d/16)*V6 + (
        231.0/16*cosine6-315.0/16*cosine4+105.0/16*cosine2-5.0/16)*V6d + (
        6435.0*(cosine4d*cosine4+cosine4*cosine4d)/128-3003.0*cosine6d/32+
        3465.0*cosine4d/64-315.0*cosine2d/32)*V8 + (6435.0/128*cosine4*cosine4
        -3003.0/32*cosine6+3465.0/64*cosine4-315.0/32*cosine2+35.0/128)*V8d + 
        (46189.0*cosine10d/256-109395.0*(cosine4d*cosine4+cosine4*cosine4d)/
        256+45045.0*cosine6d/128-15015.0*cosine4d/128+3465.0*cosine2d/256)*V10
        + (46189.0/256*cosine10-109395.0/256*cosine4*cosine4+45045.0/128*
        cosine6-15015.0/128*cosine4+3465.0/256*cosine2-63.0/256)*V10d;
    result = V0 + (3.0/2*cosine2-1.0/2)*V2 + (35.0/8*cosine4-15.0/4*cosine2+
        3.0/8)*V4 + (231.0/16*cosine6-315.0/16*cosine4+105.0/16*cosine2-5.0/16
        )*V6 + (6435.0/128*cosine4*cosine4-3003.0/32*cosine6+3465.0/64*cosine4
        -315.0/32*cosine2+35.0/128)*V8 + (46189.0/256*cosine10-109395.0/256*
        cosine4*cosine4+45045.0/128*cosine6-15015.0/128*cosine4+3465.0/256*
        cosine2-63.0/256)*V10;
}
