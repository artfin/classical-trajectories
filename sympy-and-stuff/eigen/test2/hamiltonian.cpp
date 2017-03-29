#include <chrono>

#include <iostream>
#include <math.h>
#include <Eigen/Dense>
// #include <Eigen/Core>

using namespace Eigen;
using namespace std;

const double mu1 = 14000;
const double mu2 = 38000;
const double l = 4.98; 

void inertia_tensor_dR(Matrix3d &i_dr, double &R, double &theta)
{
    i_dr(0, 0) = 2 * mu2 * R;
    i_dr(1, 0) = i_dr(2, 0) = 0;

    i_dr(0, 2) = i_dr(1, 2) = i_dr(2, 2) = 0;
    
    i_dr(0, 1) = i_dr(2, 1) = 0;
    i_dr(1, 1) = i_dr(0, 0);
}

void inertia_tensor_dtheta(Matrix3d &i_dt, double &R, double &theta)
{
    double l2 = l * l;
    double sin_t = sin(theta);
    double cos_t = cos(theta);

    i_dt(0, 0) = - mu1 * l2 * sin_t * cos_t;
    i_dt(1, 0) = 0;
    i_dt(2, 0) = - mu1 * l2 * (cos_t * cos_t - sin_t * sin_t);

    i_dt(0, 2) = i_dt(2, 0);
    i_dt(1, 2) = 0;
    i_dt(2, 2) = 2 * mu1 * l2 * sin_t * cos_t;

    i_dt(0, 1) = i_dt(2, 1) = 0;
    i_dt(1, 1) = i_dt(0, 0) + i_dt(2, 2);
}

void inertia_tensor(Matrix3d &inertia_tensor, double &R, double &theta)
{
    double sin_t = sin(theta);
    double cos_t = cos(theta);
    double l2 = l * l;
    double R2 = R * R;
    
    inertia_tensor(0, 0) = mu1 * l2 * cos_t * cos_t + mu2 * R2; 
    inertia_tensor(1, 0) = 0;
    inertia_tensor(2, 0) = -mu1 * l2 * sin_t * cos_t;

    inertia_tensor(0, 2) = inertia_tensor(2, 0);
    inertia_tensor(1, 2) = 0;
    inertia_tensor(2, 2) = mu1 * l2 * sin_t * sin_t;

    inertia_tensor(0, 1) = 0;
    inertia_tensor(1, 1) = inertia_tensor(0, 0) + inertia_tensor(2, 2);
    inertia_tensor(2, 1) = 0;

    //cout << "Inertia tensor matrix:" << endl << inertia_tensor << endl;
}

void fill_a_matrix(Matrix2d &a, double &R, double &theta)
{
    a(0, 0) = mu2;
    a(0, 1) = 0;
    a(1, 0) = 0;
    a(1, 1) = mu1 * l * l;

    //cout << "a matrix:" << endl << a << endl;
}

void fill_A_matrix(Matrix<double, 3, 2> &A, double &R, double &theta)
{
    A(0, 0) = A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = mu1 * l * l;
    A(2, 0) = A(2, 1) = 0;

    //cout << "A matrix:" << endl << A << endl;
}

double* hamiltonian(double R, double theta, double pR, double pT, double J, double alpha, double beta)
{
    //declaring angular momentum and its derivatives
    Vector3d j_vector(J * cos(alpha) * sin(beta),
                      J * sin(alpha) * sin(beta),
                      J * cos(beta));

    Vector3d j_vector_dalpha(- J * sin(alpha) * sin(beta),
                               J * cos(alpha) * sin(beta),
                               0);
    Vector3d j_vector_dbeta(J * cos(alpha) * cos(beta),
                            J * sin(alpha) * cos(beta),
                           -J * sin(beta));
    
    Vector2d p_vector(pR, pT);

    // declaring inertia tensor and it's derivatives
    Matrix<double, 3, 3> I; 
    Matrix<double, 3, 3> I_dr;
    Matrix<double, 3, 3> I_dtheta;
    
    // declaring matrices a, A
    Matrix<double, 2, 2> a;
    Matrix<double, 3, 2> A;

    // filling matrices 
    inertia_tensor(I, R, theta); 
    inertia_tensor_dR(I_dr, R, theta);
    inertia_tensor_dtheta(I_dtheta, R, theta);
    fill_a_matrix(a, R, theta);
    fill_A_matrix(A, R, theta);

    // supplementary matrices
    Matrix<double, 3, 3> I_inv = I.inverse();
    Matrix<double, 2, 2> a_inv = a.inverse();

    Matrix<double, 3, 3> G11;
    Matrix<double, 2, 2> G22;
    Matrix<double, 3, 2> G12;
   
    Matrix<double, 3, 3> t1;
    Matrix<double, 2, 2> t2;
    
    t1 = I;
    t1.noalias() -= A * a_inv * A.transpose();
    G11 = t1.inverse();
   
    t2 = a;
    t2.noalias() -= A.transpose() * I * A;
    G22 = t2.inverse();

    G12.noalias() = - G11 * A * a.inverse();
    
    //cout << "G11: " << endl << G11 << endl;
    //cout << "G22: " << endl << G22 << endl;
    //cout << "G12: " << endl << G12 << endl;

    double ang_term = 0.5 * j_vector.transpose() * G11 * j_vector;
    double kin_term = 0.5 * p_vector.transpose() * G22 * p_vector;
    double cor_term = j_vector.transpose() * G12 * p_vector;

    //cout << "ang_term: " << ang_term << endl;
    //cout << "kin_term: " << kin_term << endl;
    //cout << "cor_term: " << cor_term << endl;

    double h = ang_term = kin_term + cor_term;

    //cout << "Hamiltonian: " << h << endl;

    Matrix<double, 3, 3> G11_dr;
    Matrix<double, 3, 3> G11_dtheta;
    Matrix<double, 2, 2> G22_dr;
    Matrix<double, 2, 2> G22_dtheta;
    Matrix<double, 3, 2> G12_dr;
    Matrix<double, 3, 2> G12_dtheta;

    G11_dr = - G11 * I_dr * G11;
    G11_dtheta = - G11 * I_dtheta * G11;

    G22_dr = G22 * A.transpose() * I_inv * I_dr * I_inv * A * G22;
    G22_dtheta = G22 * A.transpose() * I_inv * I_dtheta * I_inv * A * G22;

    G12_dr = -G11_dr * A * a_inv;
    G12_dtheta = -G11_dtheta * A * a_inv;

    double ang_term_dr = 0.5 * j_vector.transpose() * G11_dr * j_vector;
    double kin_term_dr = 0.5 * p_vector.transpose() * G22_dr * p_vector;
    double cor_term_dr = j_vector.transpose() * G12_dr * p_vector;

    double h_dr = ang_term_dr + kin_term_dr + cor_term_dr;

    //cout << "Hamiltonian_dr: " << h_dr << endl;

    double ang_term_dtheta = 0.5 * j_vector.transpose() * G11_dtheta * j_vector;
    double kin_term_dtheta = 0.5 * p_vector.transpose() * G22_dtheta * p_vector;
    double cor_term_dtheta = j_vector.transpose() * G12_dtheta * p_vector;

    double h_dtheta = ang_term_dtheta + kin_term_dtheta + cor_term_dtheta;
    //cout << "Hamiltonian_dtheta: " << h_dtheta << endl;

    Vector2d h_dp = G22 * p_vector + G12.transpose() * j_vector;
    //cout << "h_dp: " << h_dp << endl;
       
    double ang_term_dalpha = j_vector.transpose() * G11 * j_vector_dalpha;
    double cor_term_dalpha = j_vector_dalpha.transpose() * G12 * p_vector;

    double h_dalpha = ang_term_dalpha + cor_term_dalpha;

    double ang_term_dbeta = j_vector.transpose() * G11 * j_vector_dbeta;
    double cor_term_dbeta = j_vector_dbeta.transpose() * G12 * p_vector;

    double h_dbeta = ang_term_dbeta + cor_term_dbeta;

    //cout << "h_dalpha: " << h_dalpha << endl;
    //cout << "h_dbeta: " << h_dbeta << endl;
   
    double *res = new double[7];

    res[0] = h;
    res[1] = h_dr;
    res[2] = h_dtheta;
    res[3] = h_dp(0);
    res[4] = h_dp(1);
    res[5] = h_dalpha;
    res[6] = h_dbeta;
    
    return res;
}

int main()
{
    double R = 1.0;
    double theta = M_PI / 2;
    double J = 10;
    double alpha = 0.05;
    double beta = 0.55;
    double pR = -10.;
    double pT = 0.1;
    
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now(); 
    
    for (int i = 0; i < 100000; i++)
    {
        double* res = hamiltonian(R, theta, pR, pT, J, alpha, beta); 
        delete[] res;
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "elapsed time: " << elapsed_seconds.count() / 100000 * pow(10, 6) << "microseconds\n";
    
    //for (int i = 0; i < 7; i++ ) {
        //cout << res[i] << endl;
    //}


    return 0;
}
