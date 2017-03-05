import sys
sys.path.append('/Users/mac/repos/sympy_project/sympy/')  # path on mac
sys.path.append('/home/artfin/Desktop/repos/sympy-project/sympy/') # path on ubuntu
sys.path.append('/home/ubuntu/sympy-project/sympy/') # path on M4-server

from lib.main import Particle, COM, Lagrange, Hamilton, LatexOutput

from sympy import *
from sympy.physics.vector import dynamicsymbols
from sympy import oo

m, M, t = symbols('m M t')
r0 = symbols('r0')

q, r1, r2 = dynamicsymbols('q r1 r2')
p, p1, p2 = dynamicsymbols('p p1 p2')

omega_x, omega_y, omega_z = symbols('omega_x omega_y omega_z')
j_x, j_y, j_z = symbols('J_x J_y J_z')

angular_velocity = [omega_x, omega_y, omega_z]
angular_momentum = [j_x, j_y, j_z]

# particle1 = Particle(m = m, x = r1 * cos(q * Rational(1, 2)), y = 0, z = - r1 * sin(q * Rational(1, 2)))
# particle2 = Particle(m = m, x = r2 * cos(q * Rational(1, 2)), y = 0, z = r2 * sin(q * Rational(1, 2)))
# particle3 = Particle(m = M, x = 0, y = 0, z = 0)

particle1 = Particle(m = m, x = -r1 * cos(q * Rational(1, 2)), y = 0, z = - r1 * sin(q * Rational(1, 2)))
particle2 = Particle(m = m, x = -r2 * cos(q * Rational(1, 2)), y = 0, z = r2 * sin(q * Rational(1, 2)))

# particles = [particle1, particle2, particle3]
particles = [particle1, particle2]

mass = sum([particle.m for particle in particles])
x_com = sum([particle.x * particle.m for particle in particles]) / mass
y_com = sum([particle.y * particle.m for particle in particles]) / mass
z_com = sum([particle.z * particle.m for particle in particles]) / mass

# for particle in particles:
# 	particle.x = particle.x - x_com
# 	particle.y = particle.y - y_com
# 	particle.z = particle.z - z_com
# 	print particle


freedom_degrees = [q, r1, r2]
freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
conjugate_momentum = [p, p1, p2]

lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives,
	angular_velocity = angular_velocity)

lagrangian = lagrange.lagrangian

eq1 = lagrangian.diff(diff(q, t)) - p
eq2 = lagrangian.diff(diff(r1, t)) - p1
eq3 = lagrangian.diff(diff(r2, t)) - p2
eq4 = lagrangian.diff(omega_x) - j_x
eq5 = lagrangian.diff(omega_y) - j_y
eq6 = lagrangian.diff(omega_z) - j_z
eqs = [eq1, eq2, eq3, eq4, eq5, eq6]

sol = solve(eqs, diff(q, t), diff(r1, t), diff(r2, t), omega_x, omega_y, omega_z)

sol[omega_x] = collect(sol[omega_x], j_x)
sol[omega_y] = collect(sol[omega_y], j_y)
sol[omega_z] = collect(sol[omega_z], j_z)

print '='*30 + '\n\n'

for key, value in sol.iteritems():
	print key, value

print '='*30 + '\n\n'

hamiltonian = lagrangian.subs([(omega_x, sol[omega_x]), (omega_y, sol[omega_y]), (omega_z, sol[omega_z]),
							   (diff(q,t), p), (diff(r1,t), p1), (diff(r2, t), p2)])
print hamiltonian



# hamilton = Hamilton(lagrange = lagrange, freedom_degrees = freedom_degrees, conjugate_momentum = conjugate_momentum)
# latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'GeneralizedWater')
