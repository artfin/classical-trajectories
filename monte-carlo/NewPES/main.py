import sys
sys.path.append("/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/NewPES/derivative_wrapper")
sys.path.append("/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/NewPES/potential_wrapper")

from potential_wrapper import potential
from derivatives_wrapper import derivative_R, derivative_Theta

import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(3, 10, 100)
u = [potential(_x, np.pi/2) for _x in x]
du_dr = [derivative_R(_x, np.pi/2) for _x in x]

axes = plt.gca()
axes.set_ylim([-1000, 1000])

plt.plot(x, u)
plt.grid()
plt.savefig('ab_initio_potential.png')


