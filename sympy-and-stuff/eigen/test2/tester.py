from ham import __rhs__
from ham import __derivatives__
import numpy as np
from time import time
import scipy.integrate as spi
import matplotlib.pyplot as plt

def extract_column(data, column):
    res = []
    for row in data:
        res.append(row[column])
    return res

temp = np.array([0,0,0,0,0,0], dtype = 'double')
def rhs(t, x):
    __rhs__(temp, *x)
    return temp

initial_conditions = np.array([10.0, np.pi/2, -1.0, 0, 1.1, 1.3])

t_start = 0.
t_end = 1e6
t_step = 10.

ode = spi.ode(rhs)

ode.set_integrator('lsoda', nsteps = 1000, method = 'bdf', atol = 1e-6)
ode.set_initial_value(initial_conditions)

start = time()

sol = []
t = [] 
while ode.successful() and ode.t < t_end:
    ode.integrate(ode.t + t_step)
    t.append(ode.t)
    sol.append(ode.y)

print 'Time needed: {0}'.format(time() - start)

R = extract_column(sol, 0)
theta = extract_column(sol, 1)

plt.plot(t, R) 
#plt.plot(t, theta)
plt.show()


