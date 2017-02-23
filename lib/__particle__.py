from sympy import *
from sympy.core.numbers import Float, Integer

class __particle__(object):
	def __init__(self, m, x, y, z):
		self.m = m
		self.x = x
		self.y = y
		self.z = z

		self.__vec__ = [self.x, self.y, self.z]

	@property
	def x(self):
		return self._x

	@property
	def y(self):
		return self._y

	@property
	def z(self):
		return self._z

	@x.setter
	def x(self, value):
		self._x = sympify(value) if isinstance(value, (int, float)) else value

	@y.setter
	def y(self, value):
		self._y = sympify(value) if isinstance(value, (int, float)) else value

	@z.setter
	def z(self, value):
		self._z = sympify(value) if isinstance(value, (int, float)) else value

	@property
	def __nvec__(self):
		_classes = (int, float, Float, Integer)
		msg = ' component has not been evulated down to int/float.'
		try:
			if not isinstance(self._nx, _classes): 
				print 'X ' + msg
				return
			if not isinstance(self._ny, _classes):
				print 'Y ' + msg 
				return
			if not isinstance(self._nz, _classes): 
				print 'Z ' + msg 
				return
		   	return [self._nx, self._ny, self._nz]
		except AttributeError:
			return None

	def __repr__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

	def __str__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

	def clear_numerical_vector(self):
		self._nx, self._ny, self._nz = None, None, None

	def evaluate(self, _subs):
		print [value for key, value in _subs.iteritems()]
		_subs = {key: value for key, value in _subs.iteritems()}

		try:
			if self._nx is not None and self._ny is not None and self._nz is not None:
				self._nx, self._ny, self._nz = [element.subs(_subs) for element in [self._nx, self._ny, self._nz]]
			else:
				self._nx, self._ny, self._nz = [element.subs(_subs) for element in self.__vec__]
		except AttributeError:
			self._nx, self._ny, self._nz = [element.subs(_subs) for element in self.__vec__]
		
	def diff(self, var, subs = None):
		if not subs:
			return [diff(self.x, var), diff(self.y, var), diff(self.z, var)]
		else:
			return [diff(self.x, var).evalf(subs = subs), 
					diff(self.y, var).evalf(subs = subs),
					diff(self.z, var).evalf(subs = subs)]