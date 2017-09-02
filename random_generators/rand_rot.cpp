#include <iostream>
#include <random>

#include <Eigen/Dense>

using namespace std;

random_device rd;
mt19937 eng( rd() );
uniform_real_distribution<double> distr( 0.0, 1.0);

void RMatrix( Matrix3d &m )
{
    double x1 = distr( eng );
    cout << "x1: " << x1 << endl;
}

int main()
{
    Matrix3d m;
    RMatrix( m );

    return 0;
}
