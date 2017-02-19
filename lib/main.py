from sympy import *
from sympy import sin, cos, diff, simplify, trigsimp, Matrix
from sympy.vector import CoordSysCartesian
from sympy.physics.vector import dynamicsymbols
from sympy import oo
from itertools import product
from pprint import pprint

from pylatex import Document, Section, Subsection, Command, Package
from pylatex.utils import italic, NoEscape

class Particle(object):
	def __init__(self, m, x, y, z):
		self.m = m
		self.x = x
		self.y = y
		self.z = z

	@property
	def vector(self):
		N = CoordSysCartesian('N')
		return self.x * N.i + self.y * N.j + self.z * N.k

	def __repr__(self):
		return 'm: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

	def __str__(self):
		return 'm: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

class Lagrange(object):
	def __init__(self, particles, freedom_degrees, freedom_degrees_derivatives, angular_velocity = None,
			inertia_transform = None, kinetic_transform = None, coriolis_transform = None):
		self.particles = particles
		self.freedom_degrees = freedom_degrees
		self.freedom_degrees_derivatives = freedom_degrees_derivatives
		self.angular_velocity = angular_velocity if angular_velocity is not None else [Symbol('omega_x'), Symbol('omega_y'), Symbol('omega_z')]

		self.inertia_transform = inertia_transform
		self.kinetic_transform = kinetic_transform
		self.coriolis_transform = coriolis_transform

		self.inertia_tensor = self.calculate_inertia_tensor()
		self.a_matrix = self.calculate_a_matrix()
		self.A_matrix = self.calculate_A_matrix()

		self.lagrangian = self.create_lagrangian()
		print 'Lagrangian: {0}'.format(self.lagrangian)

	def calculate_inertia_tensor(self):
		Ixx = [particle.m * (particle.y ** 2 + particle.z ** 2) for particle in self.particles]
		Ixx = sum(self.without_nans(Ixx))

		Iyy = [particle.m * (particle.x ** 2 + particle.z ** 2) for particle in self.particles]
		Iyy = sum(self.without_nans(Iyy))

		Izz = [particle.m * (particle.x ** 2 + particle.y ** 2) for particle in self.particles]
		Izz = sum(self.without_nans(Izz))

		Ixy = [particle.m * particle.x * particle.y for particle in self.particles]
		Ixy = - sum(self.without_nans(Ixy))

		Ixz = [particle.m * particle.x * particle.z for particle in self.particles]
		Ixz = - sum(self.without_nans(Ixz))

		Iyz = [particle.m * particle.y * particle.z for particle in self.particles]
		Iyz = - sum(self.without_nans(Iyz))

		return Matrix([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

	@staticmethod
	def without_nans(l):
		return [element for element in l if str(element) != 'nan']

	def calculate_a_matrix(self):
		"""
		Creates a matrix of needed size filled with zeros. Then calculates element by element, 
		getting rid of nans while doing this and simplifying expressions in each element.
		"""
		a = zeros(len(self.freedom_degrees), len(self.freedom_degrees))

		for j, k in product(self.freedom_degrees, self.freedom_degrees):
			interim = []
			for particle in self.particles:
				interim.append(particle.m * self.calculate_a_element(particle, j, k))

			a[self.freedom_degrees.index(j), self.freedom_degrees.index(k)] = simplify(sum(self.without_nans(interim)))
		return a

	def calculate_A_matrix(self):
		"""
		Creates a matrix of needed size filled with zeros. Then calculates element by element, 
		getting rid of nans while doing this and simplifying expressions in each element.
		"""
		A = zeros(3, len(self.freedom_degrees))

		for j, k in product(range(3), self.freedom_degrees):
			interim = []
			for particle in self.particles:
				interim.append(particle.m * Lagrange.calculate_A_element(particle, j, k))

			A[j, self.freedom_degrees.index(k)] = simplify(sum(self.without_nans(interim)))

		return A

	@staticmethod
	def calculate_a_element(particle, coord1, coord2):
		"""
		Input: particle and two generalized coordinates
		Output: simplified dot product of two derivatives (by coord1 and by coord2)
		"""
		vector = particle.vector
		derivative1 = vector.diff(coord1)
		derivative2 = vector.diff(coord2)
		return trigsimp(derivative1.dot(derivative2))

	@classmethod
	def calculate_A_element(cls, particle, number, coord):
		"""
		Input: particle, generalized coordinate and number of component needed to be returned
		(0 corresponds to X, 1 -- to Y, 2 -- to Z)
		"""
		vector = particle.vector
		derivative = vector.diff(coord)
		cross_product = vector.cross(derivative)
		return cls.return_vector_component_by_number(cross_product, number)

	@staticmethod
	def return_vector_component_by_number(vector, number):
		N = CoordSysCartesian('N')
		if number == 0: return vector.dot(N.i)
		if number == 1: return vector.dot(N.j)
		if number == 2: return vector.dot(N.k)

	def create_lagrangian(self):
		# print 'angular velocity: {0}'.format(Matrix(self.angular_velocity).shape)
		inertia_term = Rational(1, 2) * Matrix(self.angular_velocity).transpose() * self.inertia_tensor * Matrix(self.angular_velocity)
		print 'inertia term: {0}'.format(inertia_term)
		coriolis_term = Matrix(self.angular_velocity).transpose() * self.A_matrix * Matrix(self.freedom_degrees_derivatives)
		print 'coriolis term: {0}'.format(coriolis_term)
		kinetic_term = Rational(1, 2) * Matrix(self.freedom_degrees_derivatives).transpose() * self.a_matrix * Matrix(self.freedom_degrees_derivatives)
		print 'kinetic term: {0}'.format(kinetic_term)

		inertia_term = self.inertia_transform(inertia_term[0]) if self.inertia_transform is not None else inertia_term[0]
		coriolis_term = self.coriolis_transform(coriolis_term[0]) if self.coriolis_transform is not None else coriolis_term[0]
		kinetic_term = self.kinetic_transform(kinetic_term[0]) if self.kinetic_transform is not None else kinetic_term[0]

		lagrangian = inertia_term + coriolis_term + kinetic_term
		return lagrangian

class Hamilton(object):
	def __init__(self, lagrange, angular_momentum = None, conjugate_momentum = None, 
				angular_transform = None, kinetic_transform = None, coriolis_transform = None):
		self.lagrange = lagrange
		self.angular_momentum = angular_momentum if angular_momentum is not None else [Symbol('J_x'), Symbol('J_y'), Symbol('J_z')]
		self.conjugate_momentum = conjugate_momentum

		self.angular_transform = angular_transform
		self.kinetic_transform = kinetic_transform
		self.coriolis_transform = coriolis_transform

		self.G11 = (self.lagrange.inertia_tensor - self.lagrange.A_matrix * self.lagrange.a_matrix.inv() * self.lagrange.A_matrix.transpose()).inv()
		self.G22 = (self.lagrange.a_matrix - self.lagrange.A_matrix.transpose() * self.lagrange.inertia_tensor.inv() * self.lagrange.A_matrix).inv()
		self.G12 = - self.lagrange.inertia_tensor.inv() * self.lagrange.A_matrix * self.G22
		self.G21 = - self.lagrange.a_matrix.inv() * self.lagrange.A_matrix.transpose() * self.G11

		self.hamiltonian = self.create_hamiltonian()
		print 'hamiltonian: {0}'.format(self.hamiltonian)

	def create_hamiltonian(self):
		angular_term = Rational(1, 2) * Matrix(self.angular_momentum).transpose() * self.G11 * Matrix(self.angular_momentum)
		print 'angular term: {0}'.format(angular_term)
		kinetic_term = Rational(1, 2) * Matrix(self.conjugate_momentum).transpose() * self.G22 * Matrix(self.conjugate_momentum)
		print 'kinetic term: {0}'.format(kinetic_term)
		coriolis_term = Matrix(self.angular_momentum).transpose() * self.G12 * Matrix(self.conjugate_momentum)
		print 'coriolis_term: {0}'.format(coriolis_term)

		angular_term = self.angular_transform(angular_term[0]) if self.angular_transform is not None else angular_term[0]
		coriolis_term = self.coriolis_transform(coriolis_term[0]) if self.coriolis_transform is not None else coriolis_term[0]
		kinetic_term = self.kinetic_transform(kinetic_term[0]) if self.kinetic_transform is not None else kinetic_term[0]

		return angular_term + kinetic_term + coriolis_term

class COM(object):
	def __init__(self, particles):
		self.particles = particles

		self.recalculate_to_com_frame()

	@staticmethod
	def without_nans(l):
		return [element for element in l if str(element) != 'nan']

	@property
	def M(self):
		return sum([particle.m for particle in self.particles])

	@property
	def x(self):
		expr = sum(self.without_nans([particle.x * particle.m for particle in self.particles])) / self.M
		return 0 if self.M._has(oo) else expr

	@property
	def y(self):
		expr = sum(self.without_nans([particle.y * particle.m for particle in self.particles])) / self.M
		return 0 if self.M._has(oo) else expr

	@property
	def z(self):
		expr = sum(self.without_nans([particle.z * particle.m for particle in self.particles])) / self.M
		return 0 if self.M._has(oo) else expr

	def __str__(self):
		return '--- COM ---.\nX: {0};\nY: {1};\nZ: {2};'.format(self.x, self.y, self.z)

	def recalculate_to_com_frame(self):
		for particle in self.particles:
			particle.x = particle.x - self.x
			particle.y = particle.y - self.y
			particle.z = particle.z - self.z

class LatexOutput(object):
	def __init__(self, lagrangian = None, hamiltonian = None, name = 'output'):
		self.lagrangian = lagrangian
		self.hamiltonian = hamiltonian

		self.doc = Document(name)
		self.doc.packages.append(Package("breqn"))

		self.fill_document()
		self.doc.generate_pdf(clean_tex = False)
		self.doc.generate_tex()

	def fill_document(self):
		with self.doc.create(Section('Lagrangian')):
			self.doc.append(NoEscape(r'\begin{dmath}'))
			self.doc.append(NoEscape('\mathcal{L} = ' + latex(self.lagrangian)))
			self.doc.append(NoEscape(r'\end{dmath}'))

		with self.doc.create(Section('Hamiltonian')):
			self.doc.append(NoEscape(r'\begin{dmath}'))
			self.doc.append(NoEscape('\mathcal{H} = ' + latex(self.hamiltonian)))
			self.doc.append(NoEscape(r'\end{dmath}'))

# ----------------------------------------------------------

# t = Symbol('t')
# R, theta = dynamicsymbols('R theta')
# p, p_theta = dynamicsymbols('p p_theta')
# m1 = Symbol('m1')
# m2 = Symbol('m2')
# m3 = Symbol('m3')
# r0 = Symbol('r0')

# particle1 = Particle(m = m1, x = -m2/ (2 * m1 + m2 + m3) * R - r0 / 2 * cos(theta), y = r0 / 2 * sin(theta), z = 0)
# particle2 = Particle(m = m2, x = (2 * m1 + m3)  / (2 * m1 + m2 + m3) * R, y = 0, z = 0)
# particle3 = Particle(m = m1, x = - m2 / (2 * m1 + m2 + m3) * R + r0 / 2 * cos(theta), y = - r0 / 2 * sin(theta), z = 0)
# particle4 = Particle(m = m3, x = - m2 / (2 * m1 + m2 + m3) * R, y = 0, z = 0)

# particles = [particle1, particle2, particle3, particle4]

# com = COM(particles)
# print com
# particles = com.particles

# freedom_degrees = [R, theta]
# freedom_degrees_derivatives = [diff(degree, t) for degree in freedom_degrees]
# conjugate_momentum = [p, p_theta]

# lagrange = Lagrange(particles = particles, freedom_degrees = freedom_degrees, freedom_degrees_derivatives = freedom_degrees_derivatives)
# pprint(vars(lagrange))

# hamilton = Hamilton(lagrange = lagrange, conjugate_momentum = conjugate_momentum)
# pprint(vars(hamilton))

# latex_output = LatexOutput(lagrangian = lagrange.lagrangian, hamiltonian = hamilton.hamiltonian, name = 'ArCO2')




