import numpy as np
from numpy.linalg import inv
import itertools
import operator

class NumericVectorException(Exception):
	pass

class Hamiltonian(object):
	def __init__(self, particles, freedom_degrees = None, subs = None):
		self.particles = particles
		self.freedom_degrees = freedom_degrees
		self.subs = subs

		#self.check_particles()

	def check_particles(self):
		for particle in self.particles:
			if particle.__nvec__ is None:
				raise NumericVectorException

	@staticmethod
	def without_nans(l):
		return [float(element) for element in l if str(element) != 'nan']

	def calculate_inertia_tensor(self):
		Ixx = [particle.m * (particle._ny ** 2 + particle._nz ** 2) for particle in self.particles]
		Ixx = sum(self.without_nans(Ixx))

		Iyy = [particle.m * (particle._nx ** 2 + particle._nz ** 2) for particle in self.particles]
		Iyy = sum(self.without_nans(Iyy))

		Izz = [particle.m * (particle._nx ** 2 + particle._ny ** 2) for particle in self.particles]
		Izz = sum(self.without_nans(Izz))

		Ixy = [particle.m * particle._nx * particle._ny for particle in self.particles]
		Ixy = - sum(self.without_nans(Ixy))

		Ixz = [particle.m * particle._nx * particle._nz for particle in self.particles]
		Ixz = - sum(self.without_nans(Ixz))

		Iyz = [particle.m * particle._ny * particle._nz for particle in self.particles]
		Iyz = - sum(self.without_nans(Iyz))

		return np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])

	def calculate_a_matrix(self):
		"""
		Creates a matrix of needed size filled with zeros. Then calculates element by element, 
		getting rid of nans while doing this and simplifying expressions in each element.
		"""
		a = np.zeros( (len(self.freedom_degrees), len(self.freedom_degrees)) )
		for j, k in itertools.product(self.freedom_degrees, self.freedom_degrees):
			interim = []
			for particle in self.particles:
				interim.append(particle.m * self.dot_product(particle.diff(var = j, subs = self.subs), 
												 		     particle.diff(var = k, subs = self.subs)))

			a[self.freedom_degrees.index(j), self.freedom_degrees.index(k)] = sum(self.without_nans(interim))
		return a

	@staticmethod
	def dot_product(l1, l2):
		return sum(itertools.starmap(operator.mul, zip(l1, l2))) if len(l1) == len(l2) else None

	def calculate_A_matrix(self):
		"""
		Creates a matrix of needed size filled with zeros. Then calculates element by element, 
		getting rid of nans while doing this and simplifying expressions in each element.
		"""
		A = np.zeros( (3, len(self.freedom_degrees)) )

		for j, k in itertools.product(range(3), self.freedom_degrees):
			interim = []
			for particle in self.particles:
				# interim.append(particle.m * Lagrange.calculate_A_element(particle, j, k))
				# print 'nvec: {0}'.format(particle.__nvec__)
				# print 'diff_vector: {0}'.format(particle.diff(var = k, subs = self.subs))
				cross_product = np.cross(particle.__nvec__,
										 particle.diff(var = k, subs = self.subs))
				# print 'cross_product: {0}'.format(cross_product)
				interim.append(particle.m * cross_product[j])

			A[j, self.freedom_degrees.index(k)] = sum(self.without_nans(interim))

		return A

	def calculate_G_matrices(self, inertia_tensor, a_matrix, A_matrix):
		G11 = inv(inertia_tensor - A_matrix.dot(inv(a_matrix)).dot(A_matrix.transpose()))
		G22 = inv(a_matrix - A_matrix.transpose().dot(inv(inertia_tensor)).dot(A_matrix))
		G12 = - inv(inertia_tensor).dot(A_matrix).dot(G22)
		G21 = - inv(a_matrix).dot(A_matrix.transpose()).dot(G11)
		return G11, G12, G21, G22

	def calculate_hamiltonian(self, G_matrices):
		G11, G12, G21, G22 = G_matrices[0], G_matrices[1], G_matrices[2], G_matrices[3]
		J = self.subs["J"]
		p = self.subs["p"]
		theta = self.subs["theta"]
		varphi = self.subs["varphi"]
		J_vector = np.array([J*np.sin(theta)*np.cos(varphi),
							 J*np.sin(theta)*np.sin(varphi),
							 J*np.cos(theta)])
		p_vector = np.array([p])

		angular_term = 0.5 * J_vector.dot(G11).dot(J_vector.transpose())
		kinetic_term = 0.5 * p_vector.dot(G22).dot(p_vector.transpose())
		coriolis_term = J_vector.transpose().dot(G12).dot(p_vector.transpose())
		return angular_term + kinetic_term + coriolis_term