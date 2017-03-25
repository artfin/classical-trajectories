from symengine.sympy_compat import Add
from symengine import symbols
from symengine.lib.symengine_wrapper import (DenseMatrix, function_symbol, Symbol, Integer)
from symengine import sin, cos, Lambdify

from itertools import product
from time import clock

import numpy as np

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

mu1, mu2, l = symbols('mu1 mu2 l')
q1, q2 = symbols('q1 q2')

particle1 = {'x': l * sin(q1), 'y': Integer(0), 'z': l * cos(q1), 'm': mu1}
particle2 = {'x': Integer(0), 'y': Integer(0), 'z': q2, 'm': mu2}
freedom_degrees = [q1, q2]
particles = [particle1, particle2]

lagrange = Lagrange(particles, freedom_degrees)

print lagrange.I
print lagrange.a
print lagrange.A

mu1_value = 14648.
mu2_value = 38363.
l_value = 4.39

print lagrange.A[1, 0].subs({l: 1.0, mu1: 2.0})

lagrange.I = lagrange.I.subs({mu1: mu1_value, mu2: mu2_value, l: l_value})
lagrange.a = lagrange.a.subs({mu1: mu1_value, mu2: mu2_value, l: l_value})
lagrange.A = lagrange.A.subs({mu1: mu1_value, mu2: mu2_value, l: l_value})

print lagrange.I

I = Lambdify([q1, q2], lagrange.I)
a = Lambdify([q1, q2], lagrange.a)
A = Lambdify([q1, q2], lagrange.A)

attempts = int(1e5)

t0 = clock()
for i in range(attempts):
    s = I([2., 3.])
print 'Needed: {0}s'.format((clock() - t0) / attempts)

t0 = clock()
for i in range(attempts):
    s = a([2., 3.])
print 'Needed: {0}s'.format((clock() - t0) / attempts)

t0 = clock()
for i in range(attempts):
    s = A([2., 3.])
print 'Needed: {0}s'.format((clock() - t0) / attempts)
