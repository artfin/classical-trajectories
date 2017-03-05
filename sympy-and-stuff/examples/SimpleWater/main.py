import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/') # path on mac
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/') # path on ubuntu

from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

m, r0, t = symbols('m r0 t')

q = dynamicsymbols('q')
p = dynamicsymbols('p')

omega_x, omega_y, omega_z = symbols('omega_x omega_y omega_z')
j_x, j_y, j_z = symbols('j_x j_y j_z')
angular_velocity = [omega_x, omega_y, omega_z]
angular_momentum = [j_x, j_y, j_z]

particle1 = Particle(m = m, x = -r0 * cos(q * Rational(1, 2)), y = 0, z = - r0 * sin(q * Rational(1, 2)))
particle2 = Particle(m = m, x = -r0 * cos(q * Rational(1, 2)), y = 0, z = r0 * sin(q * Rational(1, 2)))
particle3 = Particle(m = oo, x = 0, y = 0, z = 0)

particles = [particle1, particle2, particle3]
com = COM(particles)
particles = com.particles

freedom_degrees = [q]
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = [p]


lagrange = Lagrange(particles = particles, 
					freedom_degrees = freedom_degrees, 
					freedom_degrees_derivatives = freedom_degrees_derivatives,
					angular_velocity = angular_velocity,)

hamilton = Hamilton(lagrange = lagrange, 
					freedom_degrees = freedom_degrees,
					conjugate_momentum = conjugate_momentum,
					angular_momentum = angular_momentum)
# latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'SimpleWater')








