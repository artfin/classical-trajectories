import numpy as np
import scipy as sp
from scipy.integrate import odeint

from pprint import pprint

def f(y, t):
	print y, t
	return y


init = [.1]
t = np.linspace(0, 1)
sol = odeint(f, init, t)














