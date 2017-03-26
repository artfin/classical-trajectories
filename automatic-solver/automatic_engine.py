# importing a wrapper for last ab-initio potential and its derivatives
import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/derivative_wrapper')
from potential_wrapper import potential
from derivatives_wrapper import derivative_R, derivative_Theta

from symengine.sympy_compat import Add
from symengine import symbols
from symengine.lib.symengine_wrapper import (DenseMatrix, function_symbol, Symbol, Integer)
from symengine import sin, cos, Lambdify

from itertools import product
from time import clock

import numpy as np
from numpy.linalg import inv

class Lagrange(object):
    """ Lagrange class symbolically constructs matrices I, a, A and their derivatives. """

    def __init__(self, particles, freedom_degrees):
        self.particles = particles
        self.freedom_degrees = freedom_degrees
        
        self.I = self.calculate_inertia_tensor()
        self.a = self.calculate_a_matrix()
        self.A = self.calculate_A_matrix()

    def calculate_inertia_tensor(self):
        Ixx = Add(*[particle['m'] * (particle['y']**2 + particle['z']**2) for particle in self.particles]).simplify()
        Iyy = Add(*[particle['m'] * (particle['x']**2 + particle['z']**2) for particle in self.particles]).simplify()
        Izz = Add(*[particle['m'] * (particle['x']**2 + particle['y']**2) for particle in self.particles]).simplify()
        Ixy = Add(*[- particle['m'] * particle['x'] * particle['y'] for particle in self.particles]).simplify()
        Ixz = Add(*[- particle['m'] * particle['x'] * particle['z'] for particle in self.particles]).simplify()
        Iyz = Add(*[- particle['m'] * particle['y'] * particle['z'] for particle in self.particles]).simplify()
             
        return DenseMatrix([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

    def calculate_a_matrix(self):

        a = DenseMatrix(len(self.freedom_degrees), len(self.freedom_degrees), [0] * len(self.freedom_degrees)**2)
        
        particle = self.particles[0]

        for j, k in product(self.freedom_degrees, self.freedom_degrees):
            index_1 = self.freedom_degrees.index(j)
            index_2 = self.freedom_degrees.index(k)
            
            res = []
            for particle in self.particles:
                derivative1 = np.array([particle['x'].diff(j), particle['y'].diff(j), particle['z'].diff(j)])
                derivative2 = np.array([particle['x'].diff(k), particle['y'].diff(k), particle['z'].diff(k)])
                res.append(particle['m'] * np.dot(derivative1, derivative2))

            a[index_1, index_2] = Add(*res).simplify()

        return a
    
    def calculate_A_matrix(self):

        A = DenseMatrix(3, len(self.freedom_degrees), [0] * len(self.freedom_degrees) * 3)

        for j, k in product(range(3), self.freedom_degrees):
            index = self.freedom_degrees.index(k)
            
            res = []
            for particle in self.particles:
                vec = np.array([particle['x'], particle['y'], particle['z']])
                derivative = np.array([particle['x'].diff(k), particle['y'].diff(k), particle['z'].diff(k)])
                res.append(particle['m'] * np.cross(vec, derivative)[j])
            
            A[j, index] = Add(*res).simplify()

        return A            

class Hamiltonian(object):
    """ Evaluates hamiltonian and it's derivatives using hard-coded automatic differentiation formualae."""

    def __init__(self, I, a, A):
        
        self.I, self.a, self.A = I, a, A

        self.I_function = Lambdify([q1, q2], self.I)
        self.a_function = Lambdify([q1, q2], self.a)
        self.A_function = Lambdify([q1, q2], self.A)
        
        I_q1_function = Lambdify([q1, q2], self.I.diff(q1))
        I_q2_function = Lambdify([q1, q2], self.I.diff(q2))
        self.Is_derivatives = np.array([I_q1_function, I_q2_function])

        a_q1_function = Lambdify([q1, q2], self.a.diff(q1))
        a_q2_function = Lambdify([q1, q2], self.a.diff(q2))
        self.as_derivatives = np.array([a_q1_function, a_q2_function])

        A_q1_function = Lambdify([q1, q2], self.A.diff(q1))
        A_q2_function = Lambdify([q1, q2], self.A.diff(q2))
        self.As_derivatives = np.array([A_q1_function, A_q2_function])

    def H_eval(self, q, p, J, theta, varphi):
        """
            Input:
                q -- np.array() of numerical values of internal coordinates
                p -- np.array() of numerical values of conjugate internal impulses
                theta, varphi -- values of spherical angles, describing the direction of J-vector

            Output: 
                H, H_derivatives
                hamiltonian numerical value and its derivatives by all q, theta and varphi
        """
        I = np.array(self.I_function(q))
        a = np.array(self.a_function(q))
        A = np.array(self.A_function(q))

        In_derivatives = []
        for derivative in self.Is_derivatives:
            In_derivatives.append(derivative(q))
        
        An_derivatives = []
        for derivative in self.As_derivatives:
            An_derivatives.append(derivative(q))

        an_derivatives = []
        for derivative in self.as_derivatives:
            an_derivatives.append(derivative(q))
        
        G11, G11_derivatives = self.G11_eval(I, a, A, In_derivatives, an_derivatives, An_derivatives)
        G22, G22_derivatives = self.G22_eval(I, a, A, In_derivatives, an_derivatives, An_derivatives)
        G12, G12_derivatives = self.G12_eval(G11, G11_derivatives, I, a, A, an_derivatives, An_derivatives)
    
        J_vector = np.array([J * np.cos(varphi) * np.sin(theta), J * np.sin(varphi) * np.sin(theta), J * np.cos(theta)])
        J_theta_derivative = np.array([J * np.cos(varphi) * np.cos(theta), J * np.sin(varphi) * np.cos(theta), - J * np.sin(theta)])
        J_varphi_derivative = np.array([- J * np.sin(varphi) * np.sin(theta), J * np.cos(varphi) * np.sin(theta), 0])

        H_value = 0.5 * np.dot(J_vector, np.dot(G11, J_vector.T)) + 0.5 * np.dot(p, np.dot(G22, p.T)) + np.dot(J_vector, np.dot(G12, p.T)) + potential(*q)
        # maybe not a very good idea but nevertheless
        potential_derivatives = [derivative_R, derivative_Theta]
        
        H_q_derivatives = np.array([])
        for G11_derivative, G22_derivative, G12_derivative, potential_derivative in zip(G11_derivatives, G22_derivatives, G12_derivatives, potential_derivatives):
            H_q_derivative = 0.5 * np.dot(J_vector, np.dot(G11_derivative, J_vector.T)) + 0.5 * np.dot(p, np.dot(G22_derivative, p.T)) + np.dot(J_vector, np.dot(G12, p.T)) + potential_derivative(*q) 
            H_q_derivatives = np.append(H_q_derivatives, H_q_derivative)
        
        H_p_derivatives = np.dot(G22, p.T) + np.dot(G12.T, J_vector.T)

        H_theta_derivative = np.dot(J_theta_derivative, np.dot(G11, J_vector)) + np.dot(J_theta_derivative, np.dot(G12, p.T))
        H_varphi_derivative = np.dot(J_varphi_derivative, np.dot(G11, J_vector)) + np.dot(J_varphi_derivative, np.dot(G12, p.T))
        

        return H_value, H_q_derivatives, H_p_derivatives, H_theta_derivative, H_varphi_derivative

    def G11_eval(self, I, a, A, In_derivatives, an_derivatives, An_derivatives):
        
        """
            Maybe useful to implement G11_eval, G22_eval, G12_eval in Cython, because they contain simple numerical operations. For-loops provide the room to speed up code.

            Input: 
                I -- numeric form of I
                a -- numeric form of a
                A -- numeric form of A
                In_derivatives -- list of numeric derivatives of I
                an_derivatives -- list of numeric derivatives of a
                An_derivatives -- list of numeric derivatives of A

            Output:
                list: G11, derivatives of G11
        """

        a_inv = inv(a)
        At = A.transpose()

        G11 = inv(I - np.dot(A, np.dot(a_inv, At)))
        
        G11_derivatives = []
        for I_derivative, a_derivative, A_derivative in zip(In_derivatives, an_derivatives, An_derivatives):
        
            b1 = A_derivative - np.dot(A, np.dot(a_inv, a_derivative))
            b2 = I_derivative - np.dot(b1, np.dot(a_inv, At))- np.dot(A, np.dot(a_inv, A_derivative.transpose()))
            G11_derivative = - np.dot(G11, np.dot(b2, G11))
            G11_derivatives.append(G11_derivative)

        return G11, G11_derivatives
    
    def G22_eval(self, I, a, A, In_derivatives, an_derivatives, An_derivatives):
        
        I_inv = inv(I)
        At = A.T

        G22 = inv(a - np.dot(At, np.dot(I_inv, A)))
        
        G22_derivatives = []
        for I_derivative, a_derivative, A_derivative in zip(In_derivatives, an_derivatives, An_derivatives):

            b1 = A_derivative.transpose() - np.dot(At, np.dot(I_inv, I_derivative))
            b2 = a_derivative - np.dot(b1, np.dot(I_inv, A)) - np.dot(At, np.dot(I_inv, A_derivative))
            G22_derivative = - np.dot(G22, np.dot(b2, G22))
            G22_derivatives.append(G22_derivative)
        
        return G22, G22_derivatives

    def G12_eval(self, G11, G11_derivatives, I, a, A, an_derivatives, An_derivatives):
        
        a_inv = inv(a)

        G12 = - np.dot(G11, np.dot(A, a_inv))

        G12_derivatives = []
        for G11_derivative, an_derivative, An_derivative in zip(G11_derivatives, an_derivatives, An_derivatives):

            b1 = np.dot(G11, np.dot(A, a_inv))
            b2 = np.dot(G11, np.dot(An_derivative, a_inv))
            b3 = np.dot(G11, np.dot(A, np.dot(a_inv, np.dot(an_derivative, a_inv)))) 
            G12_derivative = - b1 - b2 + b3
            G12_derivatives.append(G12_derivative)
        return G12_derivatives


mu1, mu2, l = symbols('mu1 mu2 l')
q1, q2 = symbols('q1 q2')

particle1 = {'x': l * sin(q1), 'y': Integer(0), 'z': l * cos(q1), 'm': mu1}
particle2 = {'x': Integer(0), 'y': Integer(0), 'z': q2, 'm': mu2}
freedom_degrees = [q1, q2]
particles = [particle1, particle2]

lagrange = Lagrange(particles, freedom_degrees)

mu1_value = 14648.
mu2_value = 38363.
l_value = 4.39

lagrange.I = lagrange.I.subs({mu1: mu1_value, mu2: mu2_value, l: l_value})
lagrange.a = lagrange.a.subs({mu1: mu1_value, mu2: mu2_value, l: l_value})
lagrange.A = lagrange.A.subs({mu1: mu1_value, mu2: mu2_value, l: l_value})

hamiltonian = Hamiltonian(lagrange.I, lagrange.a, lagrange.A)

t0 = clock()
q = np.array([1., 2.])
p = np.array([0.1, 0.2])

for i in range(10000):
    s = hamiltonian.H_eval(q, p, 10., 0.55, 0.1)
print 'Time needed" {0}s'.format((clock() - t0) / 10000)
