import sys
sys.path.append("/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/NewPES/derivative_wrapper")

from new_pes_cython import derivative_R, derivative_Theta
import matplotlib.pyplot as plt
import numpy as np

x = np.linspace(1, 10, 100)
du_dr = [derivative_R(_x, np.pi/2) for _x in x]

plt.plot(x, du_dr)
plt.show()


