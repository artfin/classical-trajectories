#include <iostream>
#include <random>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double BOLTZCONST = 1.38064e-23;
// unified atomic mass units to kg
const double DALTON = 1.660539e-27;
// atomic length unit
const double ALU = 5.29177e-11;

// reduced mass of ar and co2 = m(ar) * m(co2) / (m(ar) + m(co2)) in kg
const double MU = 20.952 * DALTON;

// hbar
const double HBAR = 1.0545718e-34;

random_device rd;
mt19937 eng( rd() );
uniform_real_distribution<double> distr( 0.0, 1.0 );

// thread_local is not supported in clang
#if defined(__clang__)
static std::mt19937 generator;
#else
static thread_local std::mt19937 generator;
#endif

// generates random rotation matrix S
void randomSMatrix( Matrix3d &m )
{
    // pick a rotation about the pole
    double theta = 2 * M_PI * distr( eng );
    // pick a direction to deflect the pole
    double phi = 2 * M_PI * distr( eng );
    // pick the amount of pole deflection
    double z = distr( eng );
    double sz = sqrt( z );

    Vector3d v ( cos(phi) * sz, sin(phi) * sz, sqrt(1 - z) );
    Matrix3d s = 2 * v * v.transpose() - Matrix<double, 3, 3>::Identity();

    Matrix3d r;
    r(0, 0) = cos(theta);
    r(0, 1) = sin(theta);
    r(0, 2) = 0;
    r(1, 0) = -sin(theta);
    r(1, 1) = cos(theta);
    r(1, 2) = 0;
    r(2, 0) = 0;
    r(2, 1) = 0;
    r(2, 2) = 1;

    // generating result
    m = s * r;
}

// normally distributed random vector ( = \dot{\vec{r}} )
Vector3d nextGaussianVec( const double &mean, const double &sigma )
{
    normal_distribution<double> d( mean, sigma );
    return Vector3d( d(generator), d(generator), d(generator) );
} 
    
int main( int argc, char* argv[] )
{
    int n = atoi( argv[1] );
    cout << "Given value of n: " << n << endl;

    Matrix3d s;
    Vector3d rdot;
    Vector3d v;
    
    const double R = 20 * ALU;
    const double R2 = pow(R, 2);
    const double temperature = 300; // K
    double sigma = sqrt( BOLTZCONST * temperature / MU );

    /*
    cout << "Boltzmann constant: " << BOLTZCONST << endl;
    cout << "temperature: " << temperature << "K " << endl;
    cout << "Dalton unit: " << DALTON << endl;
    cout << "reduced mass MU: " << MU << endl;
    cout << "sigma: " << sigma << endl;
    cout << "R: " << R << endl;
    cout << "------------------" << endl;
    */

    double omega_x, omega_y;
    double jx, jy;
    double Rdot, pR;

    double jx_au, jy_au, pR_au;

    for ( int i = 0; i < n; i++ )
    {
        // generating random rotation matrix
        randomSMatrix( s );
        // generating random vector \dot{\vec{r}}
        rdot = nextGaussianVec( 0.0, sigma );

        // v = s * rdot;
        v = rdot;
        omega_y = v(0) / R;
        omega_x = - v(1) / R;
        Rdot = v(2);

        jx = MU * R2 * omega_x;
        jy = MU * R2 * omega_y;

        pR = MU * Rdot;

        jx_au = jx / HBAR;
        jy_au = jy / HBAR;
        pR_au = pR / HBAR * ALU;

        cout << jx_au << " " << jy_au << " " << pR_au << endl;

        // cout << "jx_au: " << jx_au << endl;
        // cout << "jy_au: " << jy_au << endl;
        // cout << "pR: " << pR << endl;
        // cout << jx << " " << jy << " " << pR << endl;
    }
    
    return 0;
}
