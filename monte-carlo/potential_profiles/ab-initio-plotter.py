import sys
sys.path.append("/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/derivative_wrapper")
sys.path.append("/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper")

from potential_wrapper import potential
from derivatives_wrapper import derivative_R, derivative_Theta

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(3, 10, 100)
u = [potential(_x, np.pi/2) for _x in x]

axes = plt.gca()
axes.set_ylim([-1e-2, 1e-2])

plt.plot(x, u)
plt.grid()
plt.show()
#plt.savefig('ab-initio-potential.png')

