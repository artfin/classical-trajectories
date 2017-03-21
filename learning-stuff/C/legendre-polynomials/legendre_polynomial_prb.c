# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

# include "legendre_polynomial.h"

int main ( );
void p_exponential_product_test ( int p, double b );
void p_integral_test ( );
void p_polynomial_coefficients_test ( );
void p_polynomial_prime_test ( );
void p_polynomial_prime2_test ( );
void p_polynomial_value_test ( );
void p_polynomial_zeros_test ( );
void p_power_product_test ( int p, int e );
void p_quadrature_rule_test ( );
void pm_polynomial_value_test ( );
void pmn_polynomial_value_test ( );
void pmns_polynomial_value_test ( );
void pn_pair_product_test ( int p );
void pn_polynomial_coefficients_test ( );
void pn_polynomial_value_test ( );

/******************************************************************************/

int main ( )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for LEGENDRE_POLYNOMIAL_PRB.

  Discussion:

    LEGENDRE_POLYNOMIAL_PRB tests the LEGENDRE_POLYNOMIAL library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 March 2016

  Author:

    John Burkardt
*/
{
  double b;
  int e;
  int p;

  timestamp ( );
  printf ( "\n" );
  printf ( "LEGENDRE_POLYNOMIAL_PRB:\n" );
  printf ( "  C version.\n" );
  printf ( "  Test the LEGENDRE_POLYNOMIAL library.\n" );

  p = 5;
  b = 0.0;
  p_exponential_product_test ( p, b );

  p = 5;
  b = 1.0;
  p_exponential_product_test ( p, b );

  p_integral_test ( );

  p_polynomial_coefficients_test ( );
  p_polynomial_prime_test ( );
  p_polynomial_prime2_test ( );
  p_polynomial_value_test ( );
  p_polynomial_zeros_test ( );

  p = 5;
  e = 0;
  p_power_product_test ( p, e );

  p = 5;
  e = 1;
  p_power_product_test ( p, e );

  p_quadrature_rule_test ( );

  pm_polynomial_value_test ( );
  pmn_polynomial_value_test ( );
  pmns_polynomial_value_test ( );

  p = 5;
  pn_pair_product_test ( p );
  pn_polynomial_coefficients_test ( );
  pn_polynomial_value_test ( );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "LEGENDRE_POLYNOMIAL_PRB:\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void p_exponential_product_test ( int p, double b )

/******************************************************************************/
/*
  Purpose:

    P_EXPONENTIAL_PRODUCT_TEST tests P_EXPONENTIAL_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int P, the maximum degree of the polynomial
    factors.

    Input, double B, the coefficient of X in the exponential factor.
*/
{
  double *table;

  printf ( "\n" );
  printf ( "P_EXPONENTIAL_PRODUCT_TEST\n" );
  printf ( "  P_EXPONENTIAL_PRODUCT_TEST computes a Legendre exponential product table.\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -1.0 <= X <= +1.0 ) exp(B*X) P(I,X) P(J,X) dx\n" );
  printf ( "\n" );
  printf ( "  where P(I,X) = Legendre polynomial of degree I.\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponential argument coefficient B = %gn", b );

  table = p_exponential_product ( p, b );

  r8mat_print ( p + 1, p + 1, table, "  Exponential product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void p_integral_test ( )

/******************************************************************************/
/*
  Purpose:

    P_INTEGRAL_TEST tests P_INTEGRAL.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    15 March 2016

  Author:

    John Burkardt
*/
{
  int n;
  double value;

  printf ( "\n" );
  printf ( "P_INTEGRAL_TEST:\n" );
  printf ( "  P_INTEGRAL returns the integral of P(n,x) over [-1,+1].\n" );
  printf ( "\n" );
  printf ( "     N        Integral\n" );
  printf ( "\n" );

  for ( n = 0; n <= 10; n++ )
  {
    value = p_integral ( n );

    printf ( "  %4d  %14.6g\n", n, value );
  }

  return;
}
/******************************************************************************/

void p_polynomial_coefficients_test ( )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_COEFFICIENTS_TEST tests P_POLYNOMIAL_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int n = 10;

  printf ( "\n" );
  printf ( "P_POLYNOMIAL_COEFFICIENTS_TEST\n" );
  printf ( "  P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).\n" );

  c = p_polynomial_coefficients ( n );

  for ( i = 0; i <= n; i++ )
  {
    printf ( "\n" );
    printf ( "  P(%d,x) =\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        printf ( "%14.6g\n", c[i+j*(n+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14.6g * x\n", c[i+j*(n+1)] );
      }
      else
      {
        printf ( "%14.6g * x^%d\n", c[i+j*(n+1)], j );
      }
    }
  }
  free ( c );

  return;
}
/******************************************************************************/

void p_polynomial_prime_test ( )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_PRIME tests P_POLYNOMIAL_PRIME.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int j;
  int m;
  int n;
  double *vp;
  double *x;

  printf ( "\n" );
  printf ( "P_POLYNOMIAL_PRIME:\n" );
  printf ( "  P_POLYNOMIAL_PRIME evaluates the derivative of the\n" );
  printf ( "  Legendre polynomial P(n,x).\n" );
  printf ( "\n" );
  printf ( "                        Computed\n" );
  printf ( "     N        X           P'(N,X)\n" );
  printf ( "\n" );

  m = 11;
  a = - 1.0;
  b = + 1.0;
  x = r8vec_linspace_new ( m, a, b );

  n = 5;
  vp = p_polynomial_prime ( m, n, x );

  for ( i = 0; i < m; i++ )
  {
    printf ( "\n" );
    for ( j = 0; j <= n; j++ )
    {
      printf ( "  %4d  %12g  %24.16g\n", j, x[i], vp[i+j*m] );
    }
  }

  free ( vp );
  free ( x );

  return;
}
/******************************************************************************/

void p_polynomial_prime2_test ( )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_PRIME2_TEST tests P_POLYNOMIAL_PRIME2.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  double a;
  double b;
  int i;
  int j;
  int m;
  int n;
  double *vpp;
  double *x;

  printf ( "\n" );
  printf ( "P_POLYNOMIAL_PRIME2_TEST:\n" );
  printf ( "  P_POLYNOMIAL_PRIME2 evaluates the second derivative of the\n" );
  printf ( "  Legendre polynomial P(n,x).\n" );
  printf ( "\n" );
  printf ( "                        Computed\n" );
  printf ( "     N        X           P\"(N,X)\n" );
  printf ( "\n" );

  m = 11;
  a = - 1.0;
  b = + 1.0;
  x = r8vec_linspace_new ( m, a, b );

  n = 5;
  vpp = p_polynomial_prime2 ( m, n, x );

  for ( i = 0; i < m; i++ )
  {
    printf ( "\n" );
    for ( j = 0; j <= n; j++ )
    {
      printf ( "  %4d  %12g  %24.16g\n", j, x[i], vpp[i+j*m] );
    }
  }

  free ( vpp );
  free ( x );

  return;
}
/******************************************************************************/

void p_polynomial_value_test ( )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_VALUE_TEST tests P_POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m = 1;
  int n;
  double *v;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "P_POLYNOMIAL_VALUE_TEST:\n" );
  printf ( "  P_POLYNOMIAL_VALUE evaluates the Legendre polynomial P(n,x).\n" );
  printf ( "\n" );
  printf ( "                        Tabulated                 Computed\n" );
  printf ( "     N        X           P(N,X)                    P(N,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    p_polynomial_values ( &n_data, &n, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = p_polynomial_value ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %12g  %24.16g  %24.16g  %8e\n", n, x, fx1, fx2, e );
  }

  return;
}
/******************************************************************************/

void p_polynomial_zeros_test ( )

/******************************************************************************/
/*
  Purpose:

    P_POLYNOMIAL_ZEROS_TEST tests P_POLYNOMIAL_ZEROS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  int degree;
  double *hz;
  char title[80];
  double *z;

  printf ( "\n" );
  printf ( "P_POLYNOMIAL_ZEROS_TEST:\n" );
  printf ( "  P_POLYNOMIAL_ZEROS computes the zeros of P(n,x)\n" );
  printf ( "  Check by calling P_POLYNOMIAL_VALUE there.\n" );

  for ( degree = 1; degree <= 5; degree++ )
  {
    z = p_polynomial_zeros ( degree );
    sprintf ( title, "  Computed zeros for P(%d,z):", degree );
    r8vec_print ( degree, z, title );

    hz = p_polynomial_value ( degree, degree, z );
    sprintf ( title, "  Evaluate P(%d,z):", degree );
    r8vec_print ( degree, hz+degree*degree, title );

    free ( hz );
    free ( z );
  }
  return;
}
/******************************************************************************/

void p_power_product_test ( int p, int e )

/******************************************************************************/
/*
  Purpose:

    P_POWER_PRODUCT_TEST tests P_POWER_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int P, the maximum degree of the polynomial
    factors.

    Input, int E, the exponent of X.
*/
{
  double *table;

  printf ( "\n" );
  printf ( "P_POWER_PRODUCT_TEST\n" );
  printf ( "  P_POWER_PRODUCT_TEST computes a Legendre power product table.\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -1.0 <= X <= +1.0 ) X^E P(I,X) P(J,X) dx\n" );
  printf ( "\n" );
  printf ( "  where P(I,X) = Legendre polynomial of degree I.\n" );

  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );
  printf ( "  Exponent of X, E = %d\n", e );

  table = p_power_product ( p, e );

  r8mat_print ( p + 1, p + 1, table, "  Power product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void p_quadrature_rule_test ( )

/******************************************************************************/
/*
  Purpose:

    P_QUADRATURE_RULE_TEST tests P_QUADRATURE_RULE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  int e;
  double *f;
  int i;
  int n;
  double q;
  double q_exact;
  double *w;
  double *x;

  printf ( "\n" );
  printf ( "P_QUADRATURE_RULE_TEST:\n" );
  printf ( "  P_QUADRATURE_RULE computes the quadrature rule\n" );
  printf ( "  associated with P(n,x)\n" );

  n = 7;
  x = ( double * ) malloc ( n * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );

  p_quadrature_rule ( n, x, w );

  r8vec2_print ( n, x, w, "      X            W" );

  printf ( "\n" );
  printf ( "  Use the quadrature rule to estimate:\n" );
  printf ( "\n" );
  printf ( "    Q = Integral ( -1 <= X <= +1.0 ) X^E dx\n" );
  printf ( "\n" );
  printf ( "   E       Q_Estimate      Q_Exact\n" );
  printf ( "\n" );

  f = ( double * ) malloc ( n * sizeof ( double ) );

  for ( e = 0; e <= 2 * n - 1; e++ )
  {
    if ( e == 0 )
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = 1.0;
      }
    }
    else
    {
      for ( i = 0; i < n; i++ )
      {
        f[i] = pow ( x[i], e );
      }
    }
    q = r8vec_dot_product ( n, w, f );
    q_exact = p_integral ( e );
    printf ( "  %2d  %14.6g  %14.6g\n", e, q, q_exact );
  }

  free ( f );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void pm_polynomial_value_test ( )

/******************************************************************************/
/*
  Purpose:

    PM_POLYNOMIAL_VALUE_TEST tests PM_POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  double *v;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "PM_POLYNOMIAL_VALUE_TEST:\n" );
  printf ( "  PM_POLYNOMIAL_VALUE evaluates the Legendre polynomial Pm(n,m,x).\n" );
  printf ( "\n" );
  printf ( "                             Tabulated                 Computed\n" );
  printf ( "     N     M        X        Pm(N,M,X)                 Pm(N,M,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    pm_polynomial_values ( &n_data, &n, &m, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = pm_polynomial_value ( 1, n, m, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %4d  %12g  %24.16g  %24.16g  %8e\n", n, m, x, fx1, fx2, e );
  }

  return;
}
/******************************************************************************/

void pmn_polynomial_value_test ( )

/******************************************************************************/
/*
  Purpose:

    PMN_POLYNOMIAL_VALUE_TEST tests PMN_POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  double *v;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "PMN_POLYNOMIAL_VALUE_TEST:\n" );
  printf ( "  PMN_POLYNOMIAL_VALUE evaluates the Legendre polynomial Pmn(n,m,x).\n" );
  printf ( "\n" );
  printf ( "                             Tabulated                 Computed\n" );
  printf ( "     N     M        X       Pmn(N,M,X)                Pmn(N,M,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    pmn_polynomial_values ( &n_data, &n, &m, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = pmn_polynomial_value ( 1, n, m, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %4d  %12g  %24.16g  %24.16g  %8e\n", n, m, x, fx1, fx2, e );
  }

  return;
}
/******************************************************************************/

void pmns_polynomial_value_test ( )

/******************************************************************************/
/*
  Purpose:

    PMNS_POLYNOMIAL_VALUE_TEST tests PMNS_POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt
*/
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m;
  int n;
  double *v;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "PMNS_POLYNOMIAL_VALUE_TEST:\n" );
  printf ( "  PMNS_POLYNOMIAL_VALUE evaluates the Legendre polynomial Pmns(n,m,x).\n" );
  printf ( "\n" );
  printf ( "                             Tabulated                 Computed\n" );
  printf ( "     N     M        X       Pmns(N,M,X)                Pmns(N,M,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    pmns_polynomial_values ( &n_data, &n, &m, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = pmns_polynomial_value ( 1, n, m, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %4d  %12g  %24.16g  %24.16g  %8e\n", n, m, x, fx1, fx2, e );
  }
  return;
}
/******************************************************************************/

void pn_pair_product_test ( int p )

/******************************************************************************/
/*
  Purpose:

    PN_PAIR_PRODUCT_TEST tests PN_PAIR_PRODUCT.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int P, the maximum degree of the polynomial
    factors.
*/
{
  double *table;

  printf ( "\n" );
  printf ( "PN_PAIR_PRODUCT_TEST\n" );
  printf ( "  PN_PAIR_PRODUCT_TEST computes a pair product table for Pn(n,x).\n" );
  printf ( "\n" );
  printf ( "  Tij = integral ( -1.0 <= X <= +1.0 ) Pn(I,X) Pn(J,X) dx\n" );
  printf ( "\n" );
  printf ( "  where Pn(I,X) = normalized Legendre polynomial of degree I.\n" );
  printf ( "\n" );
  printf ( "  Maximum degree P = %d\n", p );

  table = pn_pair_product ( p );

  r8mat_print ( p + 1, p + 1, table, "  Pair product table:" );

  free ( table );

  return;
}
/******************************************************************************/

void pn_polynomial_coefficients_test ( )

/******************************************************************************/
/*
  Purpose:

    PN_POLYNOMIAL_COEFFICIENTS_TEST tests PN_POLYNOMIAL_COEFFICIENTS.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 October 2014

  Author:

    John Burkardt
*/
{
  double *c;
  int i;
  int j;
  int n = 10;

  printf ( "\n" );
  printf ( "PN_POLYNOMIAL_COEFFICIENTS_TEST\n" );
  printf ( "  PN_POLYNOMIAL_COEFFICIENTS: coefficients of normalized Legendre polynomial Pn(n,x).\n" );

  c = pn_polynomial_coefficients ( n );

  for ( i = 0; i <= n; i++ )
  {
    printf ( "\n" );
    printf ( "  P(%d,x) =\n", i );
    printf ( "\n" );
    for ( j = i; 0 <= j; j-- )
    {
      if ( c[i+j*(n+1)] == 0.0 )
      {
      }
      else if ( j == 0 )
      {
        printf ( "%14.6g\n", c[i+j*(n+1)] );
      }
      else if ( j == 1 )
      {
        printf ( "%14.6g * x\n", c[i+j*(n+1)] );
      }
      else
      {
        printf ( "%14.6g * x^%d\n", c[i+j*(n+1)], j );
      }
    }
  }
  free ( c );

  return;
}
/******************************************************************************/

void pn_polynomial_value_test ( )

/******************************************************************************/
/*
  Purpose:

    PN_POLYNOMIAL_VALUE_TEST tests PN_POLYNOMIAL_VALUE.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    18 March 2016

  Author:

    John Burkardt
*/
{
  int n_data;
  double e;
  double fx1;
  double fx2;
  double *fx2_vec;
  int m = 1;
  int n;
  double *v;
  double x;
  double x_vec[1];

  printf ( "\n" );
  printf ( "PN_POLYNOMIAL_VALUE_TEST:\n" );
  printf ( "  PN_POLYNOMIAL_VALUE evaluates the normalized Legendre polynomial Pn(n,x).\n" );
  printf ( "\n" );
  printf ( "                        Tabulated                 Computed\n" );
  printf ( "     N        X          Pn(N,X)                   Pn(N,X)                     Error\n" );
  printf ( "\n" );

  n_data = 0;

  for ( ; ; )
  {
    pn_polynomial_values ( &n_data, &n, &x, &fx1 );

    if ( n_data == 0 )
    {
      break;
    }

    x_vec[0] = x;
    fx2_vec = pn_polynomial_value ( 1, n, x_vec );
    fx2 = fx2_vec[n];
    free ( fx2_vec );

    e = fx1 - fx2;

    printf ( "  %4d  %12g  %24.16g  %24.16g  %8e\n", n, x, fx1, fx2, e );
  }

  return;
}

