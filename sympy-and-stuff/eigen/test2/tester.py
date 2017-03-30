from ham import __ham__ as hamiltonian
import numpy as np

out = np.array([0,0,0,0,0,0,0], dtype = 'double')

hamiltonian(out, 1.0, 0.55, 10.0, 2.0, 10.0, 0.55, 0.1)

print out
