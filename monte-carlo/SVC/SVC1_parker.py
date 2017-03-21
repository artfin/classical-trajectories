import autograd.numpy as np
from autograd import grad, grad_named
import vegas
from time import time
from functools import partial
import matplotlib.pyplot as plt


def parker_potential(R, theta):
    
    V0= 22.4247 * np.exp(-0.716288 * R - 0.0869136 * R**2)
    if ( R >= 6.32925 ):
        V0 -= 114.5 / R**6 + 2380. / R**8
    else:
         V0 += -0.599056 * np.exp(-0.650479 * R - 0.0320299 * R**2)
   
    V2 = 63.5744 * np.exp(-0.811806 * R - 0.075313 * R**2)
    if ( R >= 6.90026 ):
        V2 -= 26.6 / R**6 + 2080. / R**8
    else:
        V2 -= 0.207391 * np.exp(-0.620877 * R - 0.0310717 * R**2)
                                                
    V4 = 99.1128 * np.exp(-1.17577 * R -0.0477138 * R**2)
    if ( R >= 6.96450 ):
        V4  -= 410. / R**8
    else:
        V4 -= 0.0523497 * np.exp(-0.73534 * R - 0.0296750 * R**2)
    
    V6 = 318.652 * np.exp(-1.88135 * R) - 0.0374994 * np.exp(-1.05547 * R - 0.0182219 * R**2)
    V8 = 332.826 * np.exp(-2.14596 * R) - 0.0137376 * np.exp(-1.03942 * R - 0.0494781 * R**2)
    V10 = 435.837 * np.exp(-2.44616 * R) - 0.108283 * np.exp(-2.04765 * R)
    V = [V0, V2, V4, V6, V8, V10]

    cosine = np.cos(theta)
    cosine2 = np.cos(theta)**2
    cosine4 = cosine2 * cosine2
    cosine6 = cosine4 * cosine2
    cosine8 = cosine6 * cosine2
    cosine10 = cosine8 * cosine2

    L0 = 1.
    L2 = 0.5 * (3. * cosine2 - 1.)
    L4 = (35. * cosine4 - 30. * cosine2 + 3.) / 8.
    L6 = (231. * cosine6 - 315. * cosine4 + 105. * cosine2 - 5.) / 16.
    L8 = (6435. * cosine8 - 12012. * cosine6 + 6930. * cosine4 - 1260. * cosine2 + 35.) / 128.
    L10 = (46189. * cosine10 - 109395. * cosine8 + 90090. * cosine6 - 30030. * cosine4 + 3465. * cosine2 - 63.) / 256.
   
    L = [L0, L2, L4, L6, L8, L10]
        
    potential_value = sum([v * l for v, l in zip(V, L)])
    return potential_value

parker_potential_dr = grad(parker_potential, argnum = 0)
parker_potential_dtheta = grad(parker_potential, argnum = 1)

atomic_mass_unit = 1.660539040 * 10**(-27) # kg
complex_mass = (40. * 12.) / 52. * atomic_mass_unit # kg

planck_constant = 1.054571800 * 10**(-34) # joules * s
k = 1.38064852 * 10**(-23) # J / k
htoj = 4.35974417 * 10**(-18) # hartree to Joules
avogadro = 6.022 * 10**(23)
length_unit = 5.291772 * 10**(-11) # atomic units to m

def integrand(x, Temperature):
    # x = [R, theta]
    u = parker_potential(x[0], x[1]) * htoj
    du_dr = parker_potential_dr(x[0], x[1])
    
    value = np.exp(-u / (k * Temperature)) * (du_dr)**2 * x[0]**2 * np.sin(x[1])

    return value

def initialization(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[3., 20.], [0., np.pi]])
    start = time()
    result = integ(_integrand, nitn = 20, neval = 40000)
    print 'Time needed: {0}'.format(time() - start)
    print 'First integration, result = %s Q = %.2f' % (result, result.Q)

def cycle(T):
    _integrand = partial(integrand, Temperature = T)

    integ = vegas.Integrator([[3., 30.], [0., np.pi]])
    start = time()
    result = integ(_integrand, nitn = 20, neval = 5000)
	
    print 'result = %s Q = %.2f' % (result, result.Q)

    SVC_correction = np.pi * avogadro * result.mean * length_unit / 12 / (k * T)**3 * htoj**2 * planck_constant**2 / complex_mass * 10**6 # 10**6 to convert m3/mol to cm3/mol
    print 'Temperature: %d; SVC correction: %.8f' % (T, SVC_correction)
    print 'Time needed: {0}'.format(time() - start)
    
    return SVC_correction

def save_data(temperatures, svc_corrections):
    with open('data/SVC1_t_parker.dat', mode = 'w') as out:
        for temperature, svc_correction in zip(temperatures, svc_corrections):
            out.write(str(temperature) + ' ' + str(svc_correction) + '\n')

# initialization(200)
temperatures = [100 + 10 * i for i in range(2)]
svc_corrections = [cycle(temperature) for temperature in temperatures]

save_data(temperatures, svc_corrections)


