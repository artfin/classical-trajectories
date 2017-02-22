from sympy import *

class __particle__(object):
	def __init__(self, m, x, y, z):
		self.m = m
		self.x = x
		self.y = y
		self.z = z

	def __vec__(self):
		return [self.x, self.y, self.z]

	def __repr__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

	def __str__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

	def diff(self, var, subs = None):
		if not subs:
			return [diff(self.x, var), diff(self.y, var), diff(self.z, var)]
		else:
			return [diff(self.x, var).evalf(subs = subs), 
					diff(self.y, var).evalf(subs = subs),
					diff(self.z, var).evalf(subs = subs)]