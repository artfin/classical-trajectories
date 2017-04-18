#include "f.h"

double f1(const double y[], double *params)
{

  double t1,t2,t6,t7,t8,t9,t12,t13,t14,t16,t17,t18,t24,t25,t26,t32,t33,t34,t37,t38,t41,t46;
  t1 = y[0];
t2 = t1 * t1;
t6 = 0.1e1 / t2 / t1 / params[1];
t7 = y[2];
t8 = sin(t7);
t9 = t8 * t8;
t12 = cos(t7);
t13 = t12 * t12;
t14 = params[3] * params[3];
t16 = y[5];
t17 = cos(t16);
t18 = t17 * t17;
t24 = sin(t16);
t25 = y[4];
t26 = cos(t25);
t32 = t24 * t24;
t33 = t14 * t32;
t34 = t26 * t26;
t37 = y[3];
t38 = t37 * t37;
t41 = sin(t25);
t46 = t41 * t41;
return t6 / t9 * t13 * t14 * t18 + 0.2e1 * t6 / t8 * t14 * t24 * t26 * t17 * t12 + t6 * t33 * t34 + t6 * t38 - 0.2e1 * t6 * params[3] * t24 * t41 * t37 + t6 * t33 * t46;

}

double f2(const double y[],double *params)
{
  double t2,t6,t7,t8,t9,t11,t12,t13,t23,t30,t32,t33,t34,t39,t41,t42;
  t2 = y[0]*y[0];
t6 = 0.1e1 / t2 / params[1];
t7 = y[2];
t8 = sin(t7);
t9 = t8 * t8;
t11 = 0.1e1 / t9 / t8;
t12 = cos(t7);
t13 = t12 * t12;
t23 = params[2]*params[2];
t30 = params[3]*params[3];
t32 = y[5];
t33 = cos(t32);
t34 = t33 * t33;
t39 = sin(t32);
t41 = cos(y[4]);
t42 = t39 * t41;
return -(-t6 * t11 * t13 * t12 - t6 / t8 * t12 - 0.1e1 / params[0] / t23 * t11 * t12) * t30 * t34 + t6 / t9 * t30 * t42 * t33 * t13 + t6 * t30 * t42 * t33;

}

double f3(const double y[],double * params)
{
  double t2,t8,t12,t20,t23;
  t2 = params[2]*params[2];
t8 = y[0]*y[0];
t12 = 0.1e1 / t8 / params[1];
t20 = sin(y[5]);
t23 = sin(y[4]);
return 0.2e1 * (0.1e1 / t2 / params[0] / 0.2e1 + t12 / 0.2e1) * y[3] - t12 * params[3] * t20 * t23;


}

double f4(const double y[],double*params)
{
  return  y[1] / params[1];;
}

double f5(const double y[], double *params)
{
  double t2,t3,t7,t11,t12,t13,t14,t15,t16,t17,t23,t29,t30,t32,t37,t38,t40,t41,t42,t45,t49,t50,t55,t60;
  t2 = y[5];
t3 = sin(t2);
t7 = y[0]*y[0];
t11 = 0.1e1 / t7 / params[1];
t12 = y[2];
t13 = sin(t12);
t14 = t13 * t13;
t15 = 0.1e1 / t14;
t16 = cos(t12);
t17 = t16 * t16;
t23 = params[2]*params[2];
t29 = params[3];
t30 = t29 * t29;
t32 = cos(t2);
t37 = t11 / t13;
t38 = t32 * t32;
t40 = y[4];
t41 = cos(t40);
t42 = t41 * t16;
t45 = t3 * t3;
t49 = t11 * t30;
t50 = t41 * t41;
t55 = sin(t40);
t60 = t55 * t55;
return 0.1e1 / params[3] / t3 * (-0.2e1 * (t11 * t15 * t17 / 0.2e1 + 0.1e1 / params[0] / t23 * t15 / 0.2e1) * t30 * t32 * t3 + t37 * 
       t30 * t38 * t42 - t37 * t30 * t45 * t42 + t49 * t3 * t50 * t32 - t11 * t29 * t32 * t55 * y[3] + t49 * t3 * t60 * t32);

}

double f6(const double y[], double *params)
{
  double t2,t3,t7,t11,t12,t13,t15,t16,t19,t20,t22,t23,t28;
 t2 = y[5];
t3 = sin(t2);
t7 = y[0]*y[0];
t11 = 0.1e1 / t7 / params[1];
t12 = y[2];
t13 = sin(t12);
t15 = params[3];
t16 = t15 * t15;
t19 = y[4];
t20 = sin(t19);
t22 = cos(t2);
t23 = cos(t12);
t28 = cos(t19);
return -0.1e1 / params[3] / t3 * (-t11 / t13 * t16 * t3 * t20 * t22 * t23 - t11 * t15 * t3 * t28 * y[3]);

}