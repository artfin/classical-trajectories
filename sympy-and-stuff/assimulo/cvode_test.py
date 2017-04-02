from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem
import numpy as np

def rhs(t, y):
    yd = -1.0
    return np.array([yd])

y0 = [1.0]
t0 = 1.0

mod = Explicit_Problem(rhs, y0, t0)
sim = CVode(mod)

