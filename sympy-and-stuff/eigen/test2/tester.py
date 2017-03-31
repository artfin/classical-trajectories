from ham import __rhs__
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
def return_derivatives(t, x):
    #print 't: {0}; x: {1}'.format(t, x)
    __rhs__(temp, *x)
    #print 'returning: {0}'.format(temp)
    return temp

initial_conditions = np.array([10.0, 0.1, -15.0, 10., 0.1, 0.3])

t_start = 0.
t_end = 2e6
t_step = 100.

ode = spi.ode(return_derivatives)

ode.set_integrator('lsoda', nsteps = 500, method = 'bdf', atol = 1e-6)
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

plt.plot(t, R) 
plt.show()

#attempts = int(10e4)

#start = time()
#for i in range(attempts): 
    #derivatives = return_derivatives(initial_conditions)
#print 'Time needed: {0} microseconds'.format((time() - start) / attempts * 10**6)

