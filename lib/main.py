from __future__ import division

from sympy import *
from sympy import sin, cos, simplify, trigsimp, Matrix, powsimp
from sympy.simplify.cse_main import tree_cse
from sympy.vector import CoordSysCartesian
from sympy.physics.vector import dynamicsymbols
from sympy import oo
from itertools import product, combinations
from pprint import pprint

from operator import itemgetter

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
		# print 'Lagrangian: {0}'.format(self.lagrangian)

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
		# print 'inertia term: {0}'.format(inertia_term)
		coriolis_term = Matrix(self.angular_velocity).transpose() * self.A_matrix * Matrix(self.freedom_degrees_derivatives)
		# print 'coriolis term: {0}'.format(coriolis_term)
		kinetic_term = Rational(1, 2) * Matrix(self.freedom_degrees_derivatives).transpose() * self.a_matrix * Matrix(self.freedom_degrees_derivatives)
		# print 'kinetic term: {0}'.format(kinetic_term)

		inertia_term = inertia_term[0]
		coriolis_term = coriolis_term[0]
		kinetic_term = kinetic_term[0]

		lagrangian = inertia_term + coriolis_term + kinetic_term
		return lagrangian

class Hamilton(object):
	def __init__(self, lagrange, freedom_degrees = None, angular_momentum = None, conjugate_momentum = None):

		self.lagrange = lagrange
		self.angular_momentum = angular_momentum if angular_momentum is not None else [Symbol('J_x'), Symbol('J_y'), Symbol('J_z')]
		self.freedom_degrees = freedom_degrees
		self.conjugate_momentum = conjugate_momentum
		
		self.G11 = (self.lagrange.inertia_tensor - self.lagrange.A_matrix * self.lagrange.a_matrix.inv() * self.lagrange.A_matrix.transpose()).inv()
		self.G22 = (self.lagrange.a_matrix - self.lagrange.A_matrix.transpose() * self.lagrange.inertia_tensor.inv() * self.lagrange.A_matrix).inv()
		self.G12 = - self.lagrange.inertia_tensor.inv() * self.lagrange.A_matrix * self.G22
		self.G21 = - self.lagrange.a_matrix.inv() * self.lagrange.A_matrix.transpose() * self.G11

		self.hamiltonian = self.create_hamiltonian()
		print 'hamiltonian: {0}'.format(self.hamiltonian)

	@staticmethod
	def find_common_subexpression(expr):
		"""
		maybe for future use
		"""
		# preorder_traversal function allows to travel through the expression tree
		subexprs = [e for e in preorder_traversal(expr)]

		# constructing a list of dicts with count on every element got from preorder_traversal
		# make it unique
		stats = [{'element' : e, 'counts' : subexprs.count(e)} for e in subexprs]
		stats = {d['element'] : d for d in stats}.values()

		# sorting it by counts
		stats = sorted(stats, key = itemgetter('counts'), reverse = True)
		pprint(stats)

	@staticmethod
	def check_equality(expr1, expr2):
		difference = simplify(trigsimp(trigsimp(expr1 - expr2, method = 'old')))
		return True if difference == 0 else False

	def simplify_coefficients(self, coeffs):
		coeffs = [simplify(powsimp(trigsimp(trigsimp(coeff, method = 'old')))) for coeff in coeffs]

		simplified_coeffs = []

		for coeff in coeffs:
			# current complexity of coefficient
			current_complexity = count_ops(coeff)

			possible_simple_coeffs = []
			
			for degree in self.freedom_degrees:
				# implementing apart on coefficient and calculating complexity
				expr1 = apart(coeff, degree)
				complexity1 = count_ops(expr1)

				# implementing simplify on result of apart and calculating complexity
				expr2 = simplify(expr1)
				complexity2 = count_ops(expr2)

				# implementing trigsimp on result of apart and calculating complexity
				expr3 = trigsimp(expr1)
				complexity3 = count_ops(expr3)

				possible_simple_coeffs.append({'expression' : expr1, 'complexity': complexity1})
				possible_simple_coeffs.append({'expression' : expr2, 'complexity': complexity2})
				possible_simple_coeffs.append({'expression' : expr3, 'complexity': complexity3})

			# sorting all proposed simplifications in increasing order of complexity
			possible_simple_coeffs = sorted(possible_simple_coeffs, key = itemgetter('complexity'), reverse = False)
			
			# getting the best proposed simplification
			simplified_coeffs.append(possible_simple_coeffs[0]['expression']) 

		return simplified_coeffs

	@staticmethod
	def show_coefficients(combinations, coeffs):
		print '\n' + '=' * 30 + '\n'
		for term, coeff in zip(combinations, coeffs):
			print 'simplified coeff for {0}: {1}'.format(term, coeff)
		print '\n' + '=' * 30 + '\n'

	def simplify_angular_term(self, expr):
		# expanding given angular term
		expr = expand(expr)

		# angular combinations = [Jx**2, Jy**2, Jz**2, Jx*Jy, Jx*Jz, Jy*Jz]
		angular_combinations = [self.angular_momentum[i] ** 2 for i in range(3)]
		for i,j in combinations(self.angular_momentum, 2):
			angular_combinations.append(i * j)
		 
		# extracting coefficients for all angular combinations
		coeffs = [expr.coeff(term) for term in angular_combinations]

		# simplifying coefficients
		coeffs = self.simplify_coefficients(coeffs = coeffs)

		# printing coefficients for all angular terms
		self.show_coefficients(combinations = angular_combinations, coeffs = coeffs)

		# reconstructing expression
		expr_simplified = sum([term * coeff for term, coeff in zip(angular_combinations, coeffs)])

		# test equality of reconstructed expression and initial form of angular term
		difference = self.check_equality(expr, expr_simplified)
		if not difference: print 'ERROR IN RECONSTRUCTING ANGULAR EXPRESSION!'

		return expr_simplified

	def simplify_coriolis_term(self, expr):
		# expanding given coriolis term
		expr = expand(expr)

		# coriolis combinations = a list of pairs [J * p]
		coriolis_combinations = []
		for i, j in product(self.angular_momentum, self.conjugate_momentum):
			coriolis_combinations.append(i * j)

		# extracting coefficients for all coriolis combinations
		coeffs = [expr.coeff(term) for term in coriolis_combinations]

		# simplifying coefficients
		coeffs = self.simplify_coefficients(coeffs = coeffs)

		# printing coefficients for all coriolis terms
		self.show_coefficients(combinations = coriolis_combinations, coeffs = coeffs)

		# reconstructing coriolis term
		expr_simplified = sum([term * coeff for term, coeff in zip(coriolis_combinations, coeffs)])

		# test equality of reconstructed expression and initial form of coriolis term
		difference = self.check_equality(expr, expr_simplified)
		if not difference: print 'ERROR IN RECONSTRUCTING CORIOLIS TERM!'

		return expr_simplified

	def simplify_kinetic_term(self, expr):
		# expanding given kinetic term
		expr = expand(expr)

		# kinetic combinations = a list of pairs [p * p]
		kinetic_combinations = [term** 2 for term in self.conjugate_momentum]
		for i, j in combinations(self.conjugate_momentum, 2):
			kinetic_combinations.append(i * j)

		# extracting coefficients for all kinetic combinations
		coeffs = [expr.coeff(term) for term in kinetic_combinations]

		# simplifying coefficients
		coeffs = self.simplify_coefficients(coeffs = coeffs)

		# printing coefficients for all kinetic terms
		self.show_coefficients(combinations = kinetic_combinations, coeffs = coeffs)

		# reconstructing kinetic term
		expr_simplified = sum([term * coeff for term, coeff in zip(kinetic_combinations, coeffs)])

		# test equality of reconstructed expression and initial form of kinetic term
		difference = self.check_equality(expr, expr_simplified)
		if not difference: print 'ERROR IN RECONSTRUCTING CORIOLIS TERM!'

		return expr_simplified

	def create_hamiltonian(self):
		angular_term = Rational(1, 2) * Matrix(self.angular_momentum).transpose() * self.G11 * Matrix(self.angular_momentum)
		kinetic_term = Rational(1, 2) * Matrix(self.conjugate_momentum).transpose() * self.G22 * Matrix(self.conjugate_momentum)
		coriolis_term = Matrix(self.angular_momentum).transpose() * self.G12 * Matrix(self.conjugate_momentum)
		
		angular_term_simplified = self.simplify_angular_term(expr = angular_term[0])
		coriolis_term_simplified = self.simplify_coriolis_term(expr = coriolis_term[0])
		kinetic_term_simplified = self.simplify_kinetic_term(expr = kinetic_term[0])

		return angular_term_simplified + kinetic_term_simplified + coriolis_term_simplified

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
		return 0 if self.M._has(oo) else simplify(expr)

	@property
	def y(self):
		expr = sum(self.without_nans([particle.y * particle.m for particle in self.particles])) / self.M
		return 0 if self.M._has(oo) else simplify(expr)

	@property
	def z(self):
		expr = sum(self.without_nans([particle.z * particle.m for particle in self.particles])) / self.M
		return 0 if self.M._has(oo) else simplify(expr)

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
		# with self.doc.create(Section('Lagrangian')):
		# 	self.doc.append(NoEscape(r'\begin{dmath}'))
		# 	self.doc.append(NoEscape('\mathcal{L} = ' + latex(self.lagrangian)))
		# 	self.doc.append(NoEscape(r'\end{dmath}'))

		with self.doc.create(Section('Hamiltonian')):
			self.doc.append(NoEscape(r'\begin{dmath}'))
			self.doc.append(NoEscape('\mathcal{H} = ' + latex(self.hamiltonian)))
			self.doc.append(NoEscape(r'\end{dmath}'))
