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
        self.I_derivatives = np.array([I_q1_function, I_q2_function])

        a_q1_function = Lambdify([q1, q2], self.a.diff(q1))
        a_q2_function = Lambdify([q1, q2], self.a.diff(q2))
        self.a_derivatives = np.array([a_q1_function, a_q2_function])

        A_q1_function = Lambdify([q1, q2], self.A.diff(q1))
        A_q2_function = Lambdify([q1, q2], self.A.diff(q2))
        self.A_derivatives = np.array([A_q1_function, A_q2_function])

    def G11_eval(self, q, Is_derivatives, as_derivatives, As_derivatives):
        """
            Input: 
                q_values -- np.array()? of numerical values of freedom_degrees
                matrix_derivatives -- np.arrays() of symbolic matrices (I,a, A derivatives
        """

        I = np.array(self.I_function(q))
        a = np.array(self.a_function(q))
        A = np.array(self.A_function(q))

        a_inv = inv(a)
        At = A.transpose()

        G11 = inv(I - np.dot(A, np.dot(a_inv, At)))
        
        In_derivatives = []
        for derivative in Is_derivatives:
            In_derivatives.append(derivative(q))
        
        An_derivatives = []
        for derivative in As_derivatives:
            An_derivatives.append(derivative(q))

        an_derivatives = []
        for derivative in as_derivatives:
            an_derivatives.append(derivative(q))
        
        G11_derivatives = []
        for I_derivative, a_derivative, A_derivative in zip(In_derivatives, an_derivatives, An_derivatives):
        
            b1 = A_derivative - np.dot(A, np.dot(a_inv, a_derivative))
            b2 = I_derivative - np.dot(b1, np.dot(a_inv, At))- np.dot(A, np.dot(a_inv, A_derivative.transpose()))

        return G11,- np.dot(G11, np.dot(b2, G11))

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

for i in range(10000):
    s = hamiltonian.G11_eval([1., 2.], hamiltonian.I_derivatives, hamiltonian.a_derivatives, hamiltonian.A_derivatives)
print 'Time needed" {0}s'.format((clock() - t0) / 10000)
