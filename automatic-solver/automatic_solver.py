from itertools import product, chain
import autograd.numpy as np
from autograd import grad_named, grad
from autograd.numpy.linalg import inv
import scipy.integrate as spi
import json

# np.dot() instead of np.array().dot(np.array()) !

class AutomaticSolver(object):
	def __init__(self, particles, __degrees__, potential):
		self.particles = particles
		self.__degrees__ = __degrees__
		self.potential = potential

		self.dham_dq = grad(self.hamiltonian, argnum = 0)
		self.dham_dp = grad(self.hamiltonian, argnum = 1)
		self.dham_dtheta = grad(self.hamiltonian, argnum = 2)
		self.dham_dvarphi = grad(self.hamiltonian, argnum = 3)

	def dham_djx(self, q, p, theta, varphi):
		return (1. / self.J) * np.cos(theta) * np.cos(varphi) * self.dham_dtheta(q, p, theta, varphi) - \
			   (1. / self.J) * np.sin(varphi) / np.sin(theta) * self.dham_dvarphi(q, p, theta, varphi)

	def dham_djy(self, q, p, theta, varphi):
		return (1. / self.J) * np.sin(varphi) * np.cos(theta) * self.dham_dtheta(q, p, theta, varphi) + \
			   (1. / self.J) * np.cos(varphi) / np.sin(theta) * self.dham_dvarphi(q, p, theta, varphi)

	def dham_djz(self, q, p, theta, varphi):
		return - (1. / self.J) * np.sin(theta) * self.dham_dtheta(q, p, theta, varphi) 

	def hamiltonian(self, q = None, p = None, theta = None, varphi = None, effective_potential = False):
		"""
		q -- np.array object with all degrees of freedom
		p -- np.array object with conjugate momentum
		jx, jy, jz -- np.arrays of one element
		"""

		J_vector = np.array([self.J * np.cos(varphi) * np.sin(theta),
							 self.J * np.sin(varphi) * np.sin(theta), 
							 self.J * np.cos(theta)]).reshape((3,))

		Ixx = sum([particle.m * (particle.__y__(*q)**2 + particle.__z__(*q)**2) for particle in self.particles])
		Iyy = sum([particle.m * (particle.__x__(*q)**2 + particle.__z__(*q)**2) for particle in self.particles])
		Izz = sum([particle.m * (particle.__x__(*q)**2 + particle.__y__(*q)**2) for particle in self.particles])
		Ixy = - sum([particle.m * particle.__x__(*q) * particle.__y__(*q) for particle in self.particles])
		Ixz = - sum([particle.m * particle.__x__(*q) * particle.__z__(*q) for particle in self.particles])
		Iyz = - sum([particle.m * particle.__y__(*q) * particle.__z__(*q) for particle in self.particles])
		
		inertia_tensor = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
		inertia_tensor = inertia_tensor.reshape((3, 3))

		a = np.array([]).reshape((1, 0))
		for i, j in product(range(q.shape[0]), range(q.shape[0])):
			res = np.array([sum([particle.m * (particle.__dx__[i](*q) * particle.__dx__[j](*q) + \
									   		   particle.__dy__[i](*q) * particle.__dy__[j](*q) + \
									   		   particle.__dz__[i](*q) * particle.__dz__[j](*q)) 
							for particle in self.particles])]).reshape((1,1))
			a = np.hstack((a, res))

		a = a.reshape((q.shape[0], q.shape[0]))

	 	A = np.array([]).reshape((3,0))
	 	for k in range(q.shape[0]):
	 		res = np.zeros((3,1))
	 		for particle in self.particles:
	 			vec1 = np.array([particle.__x__(*q), particle.__y__(*q), particle.__z__(*q)]).reshape((3,))
	 			vec2 = np.array([particle.__dx__[k](*q), particle.__dy__[k](*q), particle.__dz__[k](*q)]).reshape((3,))
	 			res = np.hstack((res, np.array([particle.m * np.cross(vec1, vec2)]).reshape((3, 1))))

	 		A = np.hstack((A, res.sum(axis = 1).reshape((3,1))))

	 	A = A.reshape((3, q.shape[0]))

		G11 = inv(inertia_tensor - np.dot(np.dot(A, inv(a)), A.transpose()))
		G22 = inv(a - np.dot(np.dot(A.transpose(), inv(inertia_tensor)), A))
		G12 = - np.dot(np.dot(G11, A), inv(a))
		
		if not effective_potential:
			angular_component = 0.5 * np.dot(np.dot(J_vector, G11), J_vector)
			kinetic_component = 0.5 * np.dot(np.dot(p.transpose(), G22), p) 
			coriolis_component = np.dot(np.dot(J_vector, G12), p)

			return angular_component + kinetic_component + self.potential(q)
		else:
			angular_component = 0.5 * np.dot(np.dot(J_vector, inv(inertia_tensor)), J_vector)
			return angular_component + self.potential(q)


	@staticmethod
	def unpack_array(arr):
		arr_flatten = []
		for val in arr:
			if 'array' in str(type(val)):
				try:	
					arr_flatten.extend(list(chain(*val.tolist())))
				except TypeError:
					arr_flatten.extend(list(chain(val.tolist())))
			else:
				arr_flatten.append(val)
		return arr_flatten

	def pack_array(self, arr):
		return np.array(arr[0 : self.__degrees__]), \
			   np.array(arr[self.__degrees__ : 2 * self.__degrees__]), \
			   arr[2 * self.__degrees__], \
			   arr[2 * self.__degrees__ + 1]

	def rhs(self, t, y):
		print t, y

		q, p, theta, varphi = self.pack_array(y)
		vals = [q, p, theta, varphi]

		_dham_dp = self.dham_dp(*vals)
		# print 'dham_dp: {0}'.format(_dham_dp)

		_dham_dq = self.dham_dq(*vals)
		# print 'dham_dq: {0}'.format(_dham_dq)

		_dham_djx = self.dham_djx(*vals)
		_dham_djy = self.dham_djy(*vals)
		_dham_djz = self.dham_djz(*vals)

		derivatives = [_dham_dp, 
		 			  -_dham_dq, 
		 			   _dham_djx * np.sin(varphi) - _dham_djy * np.cos(varphi),
		 			  (_dham_djx * np.cos(varphi) + _dham_djy * np.sin(varphi)) * (1 / np.tan(theta)) - _dham_djz]

		# print derivatives

		return self.unpack_array(derivatives)

	def integrate(self, initial_conditions, t_start, t_end, t_step):
		"""
		initial conditions given in packed form
		"""
		init = self.unpack_array(initial_conditions)
		print 'init: {0}'.format(init)

		ode = spi.ode(self.rhs)

		ode.set_integrator('lsoda', nsteps = 500, method = 'bdf', atol = 1e-6)
		ode.set_initial_value(init, t_start)

		sol = []
		t = []
		while ode.successful() and ode.t < t_end:
			ode.integrate(ode.t + t_step)
			t.append(ode.t)
			sol.append(ode.y)

		return sol

	@staticmethod
	def extract_column(data, column):
		res = []
		for row in data:
			res.append(row[column])
		return res

	def save_file(self, filename, **kwargs):
		data = {}
		for key, value in kwargs.iteritems():
			data.update({key: self.unpack_array(value)})

		with open(filename, 'w') as out:
			out.write(json.dumps(data, indent = 4, separators = (',', ': ')))