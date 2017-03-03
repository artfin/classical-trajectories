from itertools import product
import autograd.numpy as np
from autograd.numpy.linalg import inv

def extract_column(data, column):
	res = []
	for row in data:
		res.append(row[column])
	return res

# np.dot() instead of np.array().dot(np.array()) !

####################################################################

def hamiltonian(q = None, p = None, jx = None, jy = None, jz = None, effective_potential = False):
	"""
	q -- np.array object with all degrees of freedom
	p -- np.array object with conjugate momentum
	jx, jy, jz -- np.arrays of one element
	"""
	# try:
	# 	print 'given shapes: {0}; {1}'.format(q.shape, p.shape)
	# 	print 'given values: {0}; {1}'.format(q, p)
	# except AttributeError:
	# 	pass

	J_vector = np.array([jx, jy, jz]).reshape((3,))

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
		#print 'i, j: {0}'.format(i, j)
		res = np.array([sum([particle.m * (particle.__dx__[i](*q) * particle.__dx__[j](*q) + \
								   		   particle.__dy__[i](*q) * particle.__dy__[j](*q) + \
								   		   particle.__dz__[i](*q) * particle.__dz__[j](*q)) 
						for particle in particles])]).reshape((1,1))
		a = np.hstack((a, res))

	a = a.reshape((q.shape[0], q.shape[0]))
 	# print 'a: {0}'.format(a)
 	# print 'a shape: {0}'.format(a.shape)

 	A = np.array([]).reshape((3,0))
 	for k in range(q.shape[0]):
 		res = np.zeros((3,1))
 		for particle in particles:
 			vec1 = np.array([particle.__x__(*q), particle.__y__(*q), particle.__z__(*q)]).reshape((3,))
 			vec2 = np.array([particle.__dx__[k](*q), particle.__dy__[k](*q), particle.__dz__[k](*q)]).reshape((3,))
 			# print 'vec1: {0}'.format(vec1)
 			# print 'vec1 shape: {0}'.format(vec1.shape)
 			# print 'vec2: {0}'.format(vec2)
 			# print 'cross: {0}'.format(np.cross(vec1, vec2))
 			res = np.hstack((res, np.array([particle.m * np.cross(vec1, vec2)]).reshape((3, 1))))
 			# print 'res in loop: {0}'.format(res)

 		# print 'res sum: {0}'.format(res.sum(axis = 1))
 		# print 'res sum shape: {0}'.format(res.sum(axis = 1).shape)
 		# print 'A shape: {0}'.format(A.shape)
 		A = np.hstack((A, res.sum(axis = 1).reshape((3,1))))

 	A = A.reshape((3, q.shape[0]))
 	# print 'A: {0}'.format(A)
 	# print 'A.shape: {0}'.format(A.shape)

	G11 = inv(inertia_tensor - np.dot(np.dot(A, inv(a)), A.transpose()))
	G22 = inv(a - np.dot(np.dot(A.transpose(), inv(inertia_tensor)), A))
	G12 = - np.dot(np.dot(G11, A), inv(a))

	angular_component = 0.5 * np.dot(np.dot(J_vector, G11), J_vector)
	potential = Vm / (2 * I0 * (1 - np.cos(q))) + Vp / (2 * I0 * (1 + np.cos(q)))
	
	if not effective_potential:
		kinetic_component = 0.5 * p * G22 * p 
		coriolis_component = np.dot(J_vector, G12)
		return angular_component + kinetic_component + potential
	else:
		return angular_component + potential


######################################################################
