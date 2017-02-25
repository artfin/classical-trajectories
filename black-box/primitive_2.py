import autograd.numpy as np
from autograd.core import primitive
from autograd import grad

@primitive
def logsumexp(x):
    """Numerically stable log(sum(exp(x)))"""
    max_x = np.max(x)
    return max_x + np.log(np.sum(np.exp(x - max_x)))

def logsumexp_vjp(g, ans, vs, gvs, x):
	print 'g: {0}'.format(g)
	print 'ans: {0}'.format(ans)
	print 'vs: {0}'.format(vs)
	print 'gvs: {0}'.format(gvs)
	print 'x: {0}'.format(x)
	print 'np.full g: {0}'.format(np.full(x.shape, g))
	print 'np.full ans: {0}'.format(np.full(x.shape, ans))
	return np.full(x.shape, g) * np.exp(x - np.full(x.shape, ans))

logsumexp.defvjp(logsumexp_vjp)

def example_func(y):
    z = y**2
    lse = logsumexp(z)
    return np.sum(lse)

grad_of_example = grad(example_func)
print "Gradient: ", grad_of_example(np.array([1.5, 1.0, 2.0]))
