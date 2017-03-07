from autograd import grad, grad_named
from autograd.numpy.linalg import inv
import autograd.numpy as np
from itertools import product
from time import time

class __particle__(object):
	def __init__(self, m, __x__, __y__, __z__,  vars):
		self.m = m
		self.__x__ = __x__
		self.__y__ = __y__
		self.__z__ = __z__

		self.__dx__ = [grad_named(__x__, var) for var in vars]
		self.__dy__ = [grad_named(__y__, var) for var in vars]
		self.__dz__ = [grad_named(__z__, var) for var in vars]

r0 = 1.
particle1 = __particle__(m = 1., __x__ = lambda q: r0 * np.cos(q/2), __y__ = lambda q: 0 * q, __z__ = lambda q: r0 * np.sin(q/2),
	vars = ['q'])
particle2 = __particle__(m = 1., __x__ = lambda q: r0 * np.cos(q/2), __y__ = lambda q: 0 * q, __z__ = lambda q: - r0 * np.sin(q/2),
	vars = ['q'])

particles = [particle1, particle2]

J = 10.
q = np.array([0.1])
p = np.array([0.2])
theta0 = 0.15
varphi0 = 0.01

def hamiltonian(q = None, p = None, theta = None, varphi = None):
	J_vector = np.array([J * np.cos(varphi) * np.sin(theta),
						 J * np.sin(varphi) * np.sin(theta), 
						 J * np.cos(theta)]).reshape((3,))

	Ixx = sum([particle.m * (particle.__y__(*q)**2 + particle.__z__(*q)**2) for particle in particles])
	Iyy = sum([particle.m * (particle.__x__(*q)**2 + particle.__z__(*q)**2) for particle in particles])
	Izz = sum([particle.m * (particle.__x__(*q)**2 + particle.__y__(*q)**2) for particle in particles])
	Ixy = - sum([particle.m * particle.__x__(*q) * particle.__y__(*q) for particle in particles])
	Ixz = - sum([particle.m * particle.__x__(*q) * particle.__z__(*q) for particle in particles])
	Iyz = - sum([particle.m * particle.__y__(*q) * particle.__z__(*q) for particle in particles])
	
	inertia_tensor = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
	inertia_tensor = inertia_tensor.reshape((3, 3))

	a = np.array([]).reshape((1, 0))
	for i, j in product(range(q.shape[0]), range(q.shape[0])):
		res = np.array([sum([particle.m * (particle.__dx__[i](*q) * particle.__dx__[j](*q) + \
								   		   particle.__dy__[i](*q) * particle.__dy__[j](*q) + \
								   		   particle.__dz__[i](*q) * particle.__dz__[j](*q)) 
						for particle in particles])]).reshape((1,1))
		a = np.hstack((a, res))

	a = a.reshape((q.shape[0], q.shape[0]))

 	A = np.array([]).reshape((3,0))
 	for k in range(q.shape[0]):
 		res = np.zeros((3,1))
 		for particle in particles:
 			vec1 = np.array([particle.__x__(*q), particle.__y__(*q), particle.__z__(*q)]).reshape((3,))
 			vec2 = np.array([particle.__dx__[k](*q), particle.__dy__[k](*q), particle.__dz__[k](*q)]).reshape((3,))
 			res = np.hstack((res, np.array([particle.m * np.cross(vec1, vec2)]).reshape((3, 1))))

 		A = np.hstack((A, res.sum(axis = 1).reshape((3,1))))

 	A = A.reshape((3, q.shape[0]))

	G11 = inv(inertia_tensor - np.dot(np.dot(A, inv(a)), A.transpose()))
	G22 = inv(a - np.dot(np.dot(A.transpose(), inv(inertia_tensor)), A))
	G12 = - np.dot(np.dot(G11, A), inv(a))

	angular_component = 0.5 * np.dot(np.dot(J_vector, G11), J_vector)
	kinetic_component = 0.5 * np.dot(np.dot(p.transpose(), G22), p) 
	coriolis_component = np.dot(np.dot(J_vector, G12), p)

	return angular_component + kinetic_component + coriolis_component

dham_dq = grad(hamiltonian, argnum = 0)

needed = []
for i in range(100):
	start = time()
	deriv = dham_dq(q, p, theta0, varphi0)
	needed.append(time() - start)

print 'Time needed: {0}s'.format(sum(needed) / 100)




	