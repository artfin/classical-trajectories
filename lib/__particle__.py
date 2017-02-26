from autograd import grad

class __particle__(object):

	def __init__(self, m, __x__, __y__, __z__):
		self.m = m
		self.__x__ = __x__
		self.__y__ = __y__
		self.__z__ = __z__

		self.__dx__ = grad(__x__)
		self.__dy__ = grad(__y__)
		self.__dz__ = grad(__z__)

	def __repr__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

	def __str__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)