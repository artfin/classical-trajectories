import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def g(y, x):
    y0 = y[0] # y
    y1 = y[1] # y'
    y2 = ((3*x+2)*y1 + (6*x-8)*y0)/(3*x-1)
    return y1, y2

# Initial conditions on y, y' at x=0
init = 2.0, 3.0

# First integrate from 0 to 2
x = np.linspace(0,2,100)
sol=odeint(g, init, x)

# Then integrate from 0 to -2
plt.plot(x, sol[:,0], color='b')
x = np.linspace(0,-2,100)
sol=odeint(g, init, x)
plt.plot(x, sol[:,0], color='b')

# The analytical answer in red dots
exact_x = np.linspace(-2,2,10)
exact_y = 2*np.exp(2*exact_x)-exact_x*np.exp(-exact_x)
plt.plot(exact_x,exact_y, 'o', color='r', label='exact')
plt.legend()

plt.show()