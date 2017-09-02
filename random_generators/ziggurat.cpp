#include <iostream>
#include <cmath>
#include <random>
#include <cstdlib>

#define ZIGNOR_C 128 // number of blocks
#define ZIGNOR_R 3.442619855899 // start of the right tail
#define ZIGNOR_V 9.91256303526217e-3

// DRanU() -- returns a uniform random number, U(0, 1)
// IRanU() -- returns an unsigned 32-bit integer

using namespace std;

default_random_engine generator;
uniform_real_distribution<double> distribution( 0.0, 1.0 );

// Mersenne Twister 19937 generator
random_device rd;
mt19937 eng( rd() );
uniform_int_distribution<unsigned int> distr;

// returns uniform random number
double DRanU( void )
{
    return distribution( generator );
}

// returns unsigned 32-bit integer
unsigned int IRanU( void )
{
    return distr( eng );
}

static double DRanNormalTail( double dMin, int iNegative )
{
    double x, y;
    do
    {
        x = log( DRanU() ) / dMin;
        y = log( DRanU() );
    } while( -2 * y < x * x );
    return iNegative ? x - dMin : dMin - x;
}

// s_adZigX holds coordinates, such that each rectangle as
// same are
// s_adZigR holds s_adZigX[i + 1] / s_adZigX[i]

static double s_adZigX[ZIGNOR_C + 1], s_adZigR[ZIGNOR_C];

static void zigNorInit( int iC, double dR, double dV )
{
    int i;
    double f;

    f = exp( 0.5 * dR * dR );
    s_adZigX[0] = dV / f;
    s_adZigX[1] = dR;
    s_adZigX[iC] = 0;

    for ( i = 2; i < iC; ++i )
    {
        s_adZigX[i] = sqrt( -2 * log( dV / s_adZigX[i-1] + f ));
        f = exp( -0.5 * s_adZigX[i] * s_adZigX[i] );
    }

    for ( i = 0; i < iC; ++i )
    {
        s_adZigR[i] = s_adZigX[i + 1] / s_adZigX[i];
    }
}

double DRanNormalZig( void )
{
    unsigned int i;
    double x, u, f0, f1;

    for ( ;; )
    {
        u = 2 * DRanU() - 1;
        i = IRanU() & 0x7F; // ????

        // first try the rectangular boxes
        if ( fabs(u) < s_adZigR[i] )
        {
            return u * s_adZigX[i];
        }

        // bottom box: sample from the tail
        if ( i == 0 )
        {
            return DRanNormalTail( ZIGNOR_R, u < 0 );
        }

        // is this a sample from the wedges?
        x = u * s_adZigX[i];
        f0 = exp( -0.5 * ( s_adZigX[i] * s_adZigX[i] - x * x ));
        f1 = exp( -0.5 * ( s_adZigX[i+1] * s_adZigX[i + 1] - x * x ));

        if ( f1 + DRanU() * (f0 - f1) < 1.0 )
        {
            return x;
        }
    }
}

int main( int argc, char* argv[] )
{
    int n = atoi( argv[1] ); 
    cout << "n: " << n << endl;

    for ( int i = 0; i < n; i++ )
    {
        cout << DRanNormalZig() << endl;
    }

    return 0;
}

