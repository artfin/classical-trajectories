import sys
sys.path.append("/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers/ab-initio/potential_wrapper")
sys.path.append('/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/potential_wrappers')

from potential_wrapper import potential as abpotential
from hutson import potdat4, extpot

import numpy as np
import scipy.special as sp

#import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# initialzing Hutson potential
potdat4()

A_n1 = [.224247*10**2, .635744*10**2, .991128*10**2, .318652*10**3, .332826*10**3, .435837*10**3]
A_n2 = [-.716288*10**0, -.811806*10**0, -.117577*10, -.188135*10, -.214596*10, -.244616*10] 
A_n3 = [-.869136*10**(-1), -.753131*10**(-1), -.477138*10**(-1), 0.*10**0, 0.*10**0, 0.*10**0]
B_n1 = [.599056*10**0, .207391*10**0, .523497*10**(-1), .374994*10**(-1), .137376*10**(-1), .108283*10**0]
B_n2 = [-.650479*10**0, -.620877*10**0, -.735340*10**0, -.105547*10, -.103942*10, -.204765*10]
B_n3 = [-.320299*10**(-1), -.310717*10**(-1), -.296750*10**(-1), -.182219*10**(-1), -.494781*10**(-1), 0.*10**0]
C_6 = [114.5, 26.6, 0.]
C_8 =[2380.0, 2080.0, 410.0] 
Rn = [6.32925, 6.90026, 6.96450]

def exponent_1(R, number):
	return A_n1[number] * np.exp(A_n2[number] * R + A_n3[number] * R**2)

def exponent_2(R, number):
	return B_n1[number] * np.exp(B_n2[number] * R + B_n3[number] * R**2)

def vandervaals(R, number):
	return C_6[number] / R**6 + C_8[number] / R**8

def v(R, number):
	if number < 3:
		if R > Rn[number]: 
			return exponent_1(R, number = number) - vandervaals(R, number = number)
	
	return exponent_1(R, number = number) - exponent_2(R, number = number)

L = [sp.legendre(2 * k) for k in range(6)] 

def parpotential(R, theta):
    V = [v(R, number = n) for n in range(6)]
    return sum([v(R, number = n) * L[n](np.cos(theta)) for n in range(6)])

def hutpotential(R, theta):
    return extpot(R, np.cos(theta))

def potential_well():
    fig = plt.figure()

    plt.rc('text', usetex = True)

    for index, theta in enumerate(thetavals):
    
        parvalues = [parpotential(R, theta) for R in x]
        abvalues = [abpotential(R, theta) for R in x]
        hutvalues = [hutpotential(R, theta) for R in x]

        ax = plt.subplot(2, 2, index + 1)

        ax.set_xlabel(r'\textbf{r}, (bohrs)')
        ax.set_ylabel(r'\textbf{U}, (hartree)')
        
        if index == 0:
            ax.set_title(r'Angle \theta = 0.')
        
            ax.set_xlim([7.0, 14.0])
            ax.set_ylim([-0.55e-3, 0.001e-2])
        
        if index == 1:
            ax.set_title(r'Angle $\theta = \displaystyle\frac{\pi}{6}$.')
            
            ax.set_xlim([7.0, 14.0])
            ax.set_ylim([-0.55e-3, 0.001e-2])
        
        if index == 2:
            ax.set_title(r'Angle $\theta = \displaystyle\frac{\pi}{3}$.')
        
            ax.set_xlim([6.0, 12.0])
            ax.set_ylim([-0.65e-3, 0.001e-2])

        if index == 3:
            ax.set_title(r'Angle $\theta = \displaystyle\frac{\pi}{2}$.')
        
            ax.set_xlim([5.0, 12.0])
            ax.set_ylim([-0.13e-2, 0.001e-2])


        lw = 1.5
        l1, = ax.plot(x, parvalues, color = '0.6', linestyle = 'solid', linewidth = lw)
        l2, = ax.plot(x, abvalues, color = '0.5', linestyle = 'dashed', linewidth = lw)
        l3, = ax.plot(x, hutvalues, color = '0.4', linestyle = 'dotted', linewidth = lw)

        ax.grid()

    fig.legend((l1, l2, l3), ('Parker', 'Ab-initio', 'Hutson'), 'lower center', ncol = 3, fancybox = True, shadow = True, prop = {'size': 'large'})

    plt.tight_layout(pad = 0.1, w_pad = -2.0, h_pad = -1.5)
    plt.show()

def potential_wall():
    fig = plt.figure()

    plt.rc('text', usetex = True)

    for index, theta in enumerate(thetavals):
    
        parvalues = [parpotential(R, theta) for R in x]
        abvalues = [abpotential(R, theta) for R in x]
        hutvalues = [hutpotential(R, theta) for R in x]

        ax = plt.subplot(2, 2, index + 1)

        ax.set_xlabel(r'\textbf{r}, (bohrs)')
        ax.set_ylabel(r'\textbf{U}, (hartree)')
        
        if index == 0:
            ax.set_title(r'Angle $\theta = 0$.')
            ax.set_xlim([3.0, 5.0])
            ax.set_ylim([0.0, 10.0])
        
        if index == 1:
            ax.set_title(r'Angle $\theta = \displaystyle\frac{\pi}{6}$.')
            
            ax.set_xlim([3.0, 5.0])
            ax.set_ylim([0.0, 10.0])
        
        if index == 2:
            ax.set_title(r'Angle $\theta = \displaystyle\frac{\pi}{3}$.')
        
            ax.set_xlim([3.0, 5.0])
            ax.set_ylim([0.0, 2.0])

        if index == 3:
            ax.set_title(r'Angle $\theta = \displaystyle\frac{\pi}{2}$.')
        
            ax.set_xlim([3.0, 5.0])
            ax.set_ylim([0.0, 1.0])


        lw = 1.5
        l1, = ax.plot(x, parvalues, color = '0.6', linestyle = 'solid', linewidth = lw)
        l2, = ax.plot(x, abvalues, color = '0.5', linestyle = 'dashed', linewidth = lw)
        l3, = ax.plot(x, hutvalues, color = '0.4', linestyle = 'dotted', linewidth = lw)

        ax.grid()

    fig.legend((l1, l2, l3), ('Parker', 'Ab-initio', 'Hutson'), 'lower center', ncol = 3, fancybox = True, shadow = True, prop = {'size': 'large'})

    plt.tight_layout(pad = 0.1, w_pad = -2.0, h_pad = -1.5)
    plt.show()

    plt.show()


thetavals = [0, np.pi/6, np.pi/3, np.pi/2]
x = np.linspace(1.0, 20, 500)

potential_well()
#potential_wall()
