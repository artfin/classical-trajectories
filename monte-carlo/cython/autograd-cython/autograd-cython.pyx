from autograd import grad

def f(x):
    """Simple function."""
    return x**2

gradient = grad(f)

def df(x):
    """Gradient of the simple function."""
    return gradient(x)

