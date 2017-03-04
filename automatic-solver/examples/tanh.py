import autograd.numpy as np
from autograd import grad

def tanh(x):
	print 'given x: {0} of type {1}'.format(x, type(x))
	y = np.exp(-x)
	return (1.0 - y) / (1.0 + y)

grad_tanh = grad(tanh)
print grad_tanh(1.0)






