from __future__ import print_function

import sys
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper')
from potential_wrapper import potential

import numpy as np
import scipy.special as sp
import vegas

from time import time
from functools import partial

import matplotlib.pyplot as plt
import matplotlib.cm as cm

k = 1.38064852 * 10**(-23) # J / K
htoj = 4.35974417 * 10**(-18) # hartree to J
h = 6.626070040*10**(-34)
da = 1.66 * 10**(-27) # da to kg
avogadro = 6.022 * 10**(23)
R = 8.314
patoatm = 101325

atomic_mass_unit = 9.1093826 * 10**(-31) # kg
atomic_length_unit = 5.291772 * 10**(-11) # m

def simple_integrand(x, temperature):
    # x = [R, theta]
    R = x[0]
    theta = x[1]
    potential_value = potential(R, theta)

    if potential_value < 0:
        u_kt = potential_value * htoj / (k * temperature)
        return R**2 * np.sin(theta) * sp.gammainc(2.5, -u_kt) * np.exp(-u_kt)
    else:
        return 0

def norm_arr(s):
    return [_ / max(s) for _ in s]

def simple_constant(temp):
    _simple_integrand = partial(simple_integrand, temperature = temp)
    
    print('*'*30)
    print('Temperature: {0}'.format(temp))
    
    start = time()
    integ = vegas.Integrator([[3, 20], [0, np.pi]])
    
    nitn = 9

    cmap = cm.hot

    pn = 1

    for i in [1, 5, 10, 15, 20, 25, 30, 35, 40]:
        integral = integ(_simple_integrand, nitn = i, neval = 10**4)
        print('integral (atomic units): {0}'.format(integral))
        interim = integral.mean * atomic_length_unit**3
        print('Time needed: {0}'.format(time() - start))
        
        ax = plt.subplot(3, 3, pn)
            
        ax.set_title('Iteration {0}'.format(i))
        ax.grid()
        
        x0 = []
        x1 = []
        weights = []
        for x, wgt in integ.random():
            x0.append(x[0])
            x1.append(x[1])
            weights.append(wgt)
            
        weights = norm_arr(weights)
        ax.scatter(x0, x1, color = cmap(weights), s = 1)
        
        if pn == nitn:
            break
        else:
            pn += 1

    plt.tight_layout()
    plt.show()

    eq_const = 2 * np.pi * avogadro / (R * temperature) * interim * patoatm 
    print('eq_const: {0}'.format(eq_const))
    print('*'*30 + '\n')

    return eq_const

temperatures = range(200, 201, 1)
simple_constants = [simple_constant(temperature) for temperature in temperatures]

