from autograd import grad
from autograd import grad_named

class __particle__(object):

	def __init__(self, m, __x__, __y__, __z__,  vars):

		self.m = m
		self.__x__ = __x__
		self.__y__ = __y__
		self.__z__ = __z__

		self.__dx__ = [grad_named(__x__, var) for var in vars]
		self.__dy__ = [grad_named(__y__, var) for var in vars]
		self.__dz__ = [grad_named(__z__, var) for var in vars]

	def __repr__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

	def __str__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)