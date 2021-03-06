import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/sympy-and-stuff')  # path on mac
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/') # path on ubuntu
sys.path.append('/home/ubuntu/sympy-project/sympy/') # path on M4-server

from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

m, M, t = symbols('m M t')

q, r1, r2 = dynamicsymbols('q r1 r2')
p, p1, p2 = dynamicsymbols('p p1 p2')

particle1 = Particle(m = m, x = -r1 * cos(q * Rational(1, 2)), y = 0, z = - r1 * sin(q * Rational(1, 2)))
particle2 = Particle(m = m, x = -r2 * cos(q * Rational(1, 2)), y = 0, z = r2 * sin(q * Rational(1, 2)))
particle3 = Particle(m = M, x = 0, y = 0, z = 0)
particles = [particle1, particle2, particle3]

# list of freedom degrees
freedom_degrees = [q, r1, r2]
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = [p, p1, p2]

com = COM(particles, freedom_degrees = freedom_degrees)
print com
x_com = com.x
y_com = com.y
z_com = com.z

particle1.x = particle1.x - x_com
particle1.y = particle1.y - y_com
particle1.z = particle1.z - z_com
print 'PARTICLE 1: {0}\n\n'.format(particle1)

particle2.x = particle2.x - x_com
particle2.y = particle2.y - y_com
particle2.z = particle2.z - z_com
print 'PARTICLE 2: {0}\n'.format(particle2)

particle3.x = particle3.x - x_com
particle3.y = particle3.y - y_com
particle3.z = particle3.z - z_com
print 'PARTICLE 3: {0}\n'.format(particle3)



# lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives)
# hamilton = Hamilton(lagrange = lagrange, freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)
# latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'FullWater')
