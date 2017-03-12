from sympy import symbols, Rational, simplify, zeros, S, diff, solve, trigsimp
from sympy.physics.vector import dynamicsymbols
from sympy import sin, cos, Matrix
from pprint import pprint
from itertools import product
from time import time
import dill
import pickle

def calculate_inertia_tensor(particles):
	Ixx = simplify(sum([particle['m'] * (particle['y']**2 + particle['z'] ** 2) for particle in particles]))
	Iyy = simplify(sum([particle['m'] * (particle['x']**2 + particle['z'] ** 2) for particle in particles]))
	Izz = simplify(sum([particle['m'] * (particle['x']**2 + particle['y'] ** 2) for particle in particles]))
	Ixy = simplify(sum([particle['m'] * particle['x'] * particle['y'] for particle in particles]))
	Ixz = simplify(sum([particle['m'] * particle['x'] * particle['z'] for particle in particles]))
	Iyz = simplify(sum([particle['m'] * particle['y'] * particle['z'] for particle in particles]))

	return Matrix([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

def calculate_A_matrix(particles, freedom_degrees):
	A = zeros(3, len(freedom_degrees))

	for j, k in product(range(3), freedom_degrees):
		interim = []
		for particle in particles:
			v1 = Matrix([particle['x'], particle['y'], particle['z']])
			v2 = Matrix([particle['x'].diff(k), particle['y'].diff(k), particle['z'].diff(k)])
			cross_product = v1.cross(v2)

			interim.append(particle['m'] * cross_product[j])

		A[j, freedom_degrees.index(k)] = simplify(sum(interim))

	return A

def calculate_a_matrix(particles, freedom_degrees):
	a = zeros(len(freedom_degrees), len(freedom_degrees))

	for j, k in product(freedom_degrees, freedom_degrees):
		interim = []
		for particle in particles:
			v1 = Matrix([particle['x'].diff(j), particle['y'].diff(j), particle['z'].diff(j)])
			v2 = Matrix([particle['x'].diff(k), particle['y'].diff(k), particle['z'].diff(k)])

			interim.append(particle['m'] * v1.dot(v2))

		a[freedom_degrees.index(j), freedom_degrees.index(k)] = simplify(sum(interim))
	return a

m1, m2, m3, r0, t = symbols('m1 m2 m3 r0 t')
R, theta = dynamicsymbols('R theta')

omega_x, omega_y, omega_z = symbols('Omega_x Omega_y Omega_z')
angular_velocity = Matrix([omega_x, omega_y, omega_z])

particle1 = {'m': m1,
			 'x': -m2 / (2 * m1 + m2 + m3) * R - r0 * Rational(1, 2) * cos(theta),
			 'y': r0 * Rational(1, 2) * sin(theta),
			 'z': S(0)}

particle2 = {'m': m2,
			 'x': (2 * m1 + m3) / (2 * m1 + m2 + m3) * R,
			 'y': S(0),
			 'z': S(0)}

particle3 = {'m': m1,
			 'x': - m2 / (2 * m1 + m2 + m3) * R + r0 * Rational(1, 2) * cos(theta), 
			 'y': - r0 * Rational(1, 2) * sin(theta),
			 'z': S(0)}

particle4 = {'m': m3,
			 'x': - m2 / (2 * m1 + m2 + m3) * R,
			 'y': S(0),
			 'z': S(0)}

particles = [particle1, particle2, particle3, particle4]

freedom_degrees = [R, theta]

freedom_degrees_deriv = Matrix([diff(R, t), diff(theta, t)])

inertia_tensor = calculate_inertia_tensor(particles)
a_matrix = calculate_a_matrix(particles, freedom_degrees)
A_matrix = calculate_A_matrix(particles, freedom_degrees)

angular_component = angular_velocity.dot(inertia_tensor.dot(angular_velocity.transpose()))
coriolis_component = angular_velocity.dot(A_matrix.dot(freedom_degrees_deriv.transpose()))
kinetic_component = freedom_degrees_deriv.dot(a_matrix.dot(freedom_degrees_deriv.transpose()))
lagrangian = angular_component + coriolis_component + kinetic_component

# substituting euler angles

alpha, phi, psi = dynamicsymbols('alpha phi psi') # theta, phi, psi
V = Matrix([[sin(alpha) * sin(psi), cos(psi), 0], [sin(alpha) * cos(psi), - sin(psi), 0], [cos(alpha), 0, 1]])
e_dot = Matrix([diff(alpha, t), diff(phi, t), diff(psi, t)])

lagrangian_R = diff(lagrangian, R)
lagrangian_R_dot = diff(lagrangian, diff(R, t))
lagrangian_theta = diff(lagrangian, theta)
lagrangian_theta_dot = diff(lagrangian, diff(theta, t))

eq1 = diff(lagrangian_R_dot, t) - lagrangian_R 
eq2 = diff(lagrangian_theta_dot, t) - lagrangian_theta

lagrangian = lagrangian.subs({omega_x: V.dot(e_dot)[0],
						      omega_y: V.dot(e_dot)[1],
						      omega_z: V.dot(e_dot)[2]})

lagrangian_alpha = diff(lagrangian, alpha)
lagrangian_alpha_dot = diff(lagrangian, diff(alpha, t))
lagrangian_phi = diff(lagrangian, phi)
lagrangian_phi_dot = diff(lagrangian, diff(phi, t))
lagrangian_psi = diff(lagrangian, psi)
lagrangian_psi_dot = diff(lagrangian, diff(psi, t))

eq3 = diff(lagrangian_alpha_dot, t) - lagrangian_alpha
eq4 = diff(lagrangian_phi_dot, t) - lagrangian_phi
eq5 = diff(lagrangian_psi_dot, t) - lagrangian_psi

start = time()
R_2dot = solve(eq1, diff(R, t, t))
print 'R_2dot. Needed: {0}'.format(time() - start)
print '-'*20

start = time()
theta_2dot = solve(eq2, diff(theta, t, t))
print 'theta_2dot. Needed: {0}'.format(time() - start)
print '-'*20

start = time()
alpha_2dot = solve(eq3, diff(alpha, t, t))
print 'alpha_2dot. Needed: {0}'.format(time() - start)
print '-'*20

start = time()
phi_2dot = solve(eq4, diff(phi, t, t))
print 'phi_2dot. Needed: {0}'.format(time() - start)
print '-'*20

start = time()
psi_2dot = solve(eq5, diff(psi, t, t))
print 'psi_2dot. Needed: {0}'.format(time() - start)
print '-'*20

start = time()
theta_dot = solve([alpha_2dot, phi_2dot, psi_2dot], diff(theta, t))
print 'theta_dot. Needed: {0}'.format(time() - start)
print '-'*20

start = time()
phi_dot = solve([alpha_2dot, phi_2dot, psi_2dot], diff(phi, t))
print 'phi_dot. Needed: {0}'.format(time() - start)
print '-'*20

start = time()
psi_dot = solve([alpha_2dot, phi_2dot, psi_2dot], diff(psi, t))
print 'psi_dot. Needed: {0}'.format(time() - start)
print '-'*20

