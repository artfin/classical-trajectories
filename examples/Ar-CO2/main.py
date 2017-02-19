import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/') # path on mac
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/') # path on ubuntu

from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

m1, m2, m3, r0, t = symbols('m1 m2 m3 r0 t')

R, theta = dynamicsymbols('R theta')
p, p_theta = dynamicsymbols('p p_theta')

particle1 = Particle(m = m1, x = -m2 / (2 * m1 + m2 + m3) * R - r0 * Rational(1, 2) * cos(theta), y = r0 * Rational(1, 2) * sin(theta), z = 0)
particle2 = Particle(m = m2, x = (2 * m1 + m3)  / (2 * m1 + m2 + m3) * R, y = 0, z = 0)
particle3 = Particle(m = m1, x = - m2 / (2 * m1 + m2 + m3) * R + r0 * Rational(1, 2) * cos(theta), y = - r0 * Rational(1, 2) * sin(theta), z = 0)
particle4 = Particle(m = m3, x = - m2 / (2 * m1 + m2 + m3) * R, y = 0, z = 0)

particles = [particle1, particle2, particle3, particle4]
com = COM(particles)
particles = com.particles

freedom_degrees = [R, theta]
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = [p, p_theta]


lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives)
hamilton = Hamilton(lagrange = lagrange, conjugate_momentum = conjugate_momentum)
latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'Ar-CO2')
