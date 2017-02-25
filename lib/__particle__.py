class __particle__(object):

	def __init__(self, m, x, y, z, dx, dy, dz):
		self.m = m
		self.x = x
		self.y = y
		self.z = z

		self.dx = dx
		self.dy = dy
		self.dz = dz

	def __vec__(self, q):
		return [self.x(q), self.y(q), self.z(q)]

	def __dvec__(self, q):
		return [self.dx(q), self.dy(q), self.dz(q)]

	def __repr__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)

	def __str__(self):
		return 'mass: {0}; x: {1}; y: {2}; z: {3}'.format(self.m, self.x, self.y, self.z)