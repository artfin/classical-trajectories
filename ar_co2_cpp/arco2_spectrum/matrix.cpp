#include <iostream>
#include <math.h>
#include <Eigen/Dense>

#include "../parker_snow/psp_pes.h"
// #include <Eigen/Core>
#include "matrix.h"

using namespace Eigen;
using namespace std;

const double mu1 = 14579.0;
const double mu2 = 38183.0;
const double l = 4.398; 

void inertia_tensor_dR(Matrix3d &i_dr, double &R, double &theta)
{
    i_dr(0, 0) = i_dr(1, 1) = 2 * mu2 * R;
    
    i_dr(1, 0) = i_dr(2, 0) = 0;    
    i_dr(0, 2) = i_dr(1, 2) = i_dr(2, 2) = 0;
    i_dr(0, 1) = i_dr(2, 1) = 0;
}

void inertia_tensor_dtheta(Matrix3d &i_dt, double &R, double &theta)
{
    double l2 = l * l;
    double sin_t = sin(theta);
    double cos_t = cos(theta);

    i_dt(0, 0) = - 2 * mu1 * l2 * sin_t * cos_t;
    i_dt(1, 0) = 0;
    i_dt(2, 0) = i_dt(0, 2) = - mu1 * l2 * (cos_t * cos_t - sin_t * sin_t);

    i_dt(1, 2) = 0;
    i_dt(2, 2) = - i_dt(0, 0); 

    i_dt(0, 1) = i_dt(1, 1) = i_dt(2, 1) = 0;
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
}

void fill_a_matrix(Matrix2d &a, double &R, double &theta)
{
    a(0, 0) = mu2;
    a(0, 1) = 0;
    a(1, 0) = 0;
    a(1, 1) = mu1 * l * l;
}

void fill_A_matrix(Matrix<double, 3, 2> &A, double &R, double &theta)
{
    A(0, 0) = A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = mu1 * l * l;
    A(2, 0) = A(2, 1) = 0;
}

void hamiltonian(double* out, double R, double theta, double pR, double pT, double alpha, double beta, double J)
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
    t2.noalias() -= A.transpose() * I_inv * A;
    G22 = t2.inverse();

    G12.noalias() = - G11 * A * a.inverse();
    
    Matrix<double, 3, 3> G11_dr;
    Matrix<double, 3, 3> G11_dtheta;
    Matrix<double, 2, 2> G22_dr;
    Matrix<double, 2, 2> G22_dtheta;
    Matrix<double, 3, 2> G12_dr;
    Matrix<double, 3, 2> G12_dtheta;
    

    // derivatives dH/dr and dH/dtheta
    G11_dr = - G11 * I_dr * G11;
    G11_dtheta = - G11 * I_dtheta * G11;
   
    G22_dr = - G22 * A.transpose() * I_inv * I_dr * I_inv * A * G22; 
    G22_dtheta = - G22 * A.transpose() * I_inv * I_dtheta * I_inv * A * G22;

    G12_dr = - G11_dr * A * a_inv;
    G12_dtheta = -G11_dtheta * A * a_inv;

    double ang_term, kin_term, cor_term;

    ang_term = 0.5 * j_vector.transpose() * G11_dr * j_vector;
    kin_term = 0.5 * p_vector.transpose() * G22_dr * p_vector;
    cor_term = j_vector.transpose() * G12_dr * p_vector;
    double h_dr = ang_term + kin_term + cor_term;
    
    //cout << "h_dr: " << h_dr << endl;

    ang_term = 0.5 * j_vector.transpose() * G11_dtheta * j_vector;
    kin_term = 0.5 * p_vector.transpose() * G22_dtheta * p_vector;
    cor_term = j_vector.transpose() * G12_dtheta * p_vector;
    double h_dtheta = ang_term + kin_term + cor_term;
    
    // (vector of) derivatives dH/dp
    Vector2d h_dp = G22 * p_vector + G12.transpose() * j_vector;

    // derivatives dH/dalpha, dH/dbeta
    ang_term = j_vector.transpose() * G11 * j_vector_dalpha;
    cor_term = j_vector_dalpha.transpose() * G12 * p_vector;
    double h_dalpha = ang_term + cor_term;

    ang_term = j_vector.transpose() * G11 * j_vector_dbeta;
    cor_term = j_vector_dbeta.transpose() * G12 * p_vector;
    double h_dbeta = ang_term + cor_term;
    
    out[0] = h_dr; 
    out[1] = h_dtheta;
    out[2] = h_dp(0);
    out[3] = h_dp(1);
    out[4] = h_dalpha;
    out[5] = h_dbeta;

    //cout << "out[0]: " << out[0] << endl;
}




void rhs(double* out, double R, double theta, double pR, double pT, double alpha, double beta, double J)
// input:
//      out -- prepared array to fill the right-hand sides of ODEs
{
    double Jsint = J * sin(theta);

    double* derivatives = new double[6];
    hamiltonian(derivatives, R, theta, pR, pT, alpha, beta, J);
    
    out[0] = derivatives[2]; // /dot(R) = dH/dpR
    out[1] = derivatives[3]; // /dot(theta) = dH/dpT
    out[2] = - derivatives[0] - dpsp_pesdR(R,theta);  // /dot(pR) = - dH/dR = - dT/dR - dU/dR (to hartrees from cm^-1)
    out[3] = - derivatives[1] - dpsp_pesdTheta(R, theta); // /dot(pT) = - dH/dtheta = -dT/dtheta - dU/dtheta (to hartrees from cm^-1)
    out[4] = 1 / Jsint * derivatives[5]; // /dot(varphi) = 1 / J / sin(teta) * dH/dtheta
    out[5] = - 1 / Jsint * derivatives[4]; // /dot(theta) = - 1 / J / sin(theta) * dH/dvarphi

    delete[] derivatives;
}
