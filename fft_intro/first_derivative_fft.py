import numpy as np
import matplotlib.pyplot as plt

N = 50
xmin = 0
xmax = 2.0 * np.pi
step = (xmax - xmin) / N

xdata = np.linspace(step, xmax, N)
print('xdata: {0}'.format(xdata))

v = np.exp( np.sin( xdata ))
derv = np.cos(xdata) * v 
vhat = np.fft.fft(v)

what = 1j * np.zeros(N)
what[0 : (N / 2)] = 1j * np.arange(0, N/2, 1)
what[(N / 2) + 1:] = 1j * np.arange(-N/2 + 1.0, 0, 1)
print('what: {0}'.format(what))

what = what * vhat
w = np.real( np.fft.ifft(what) )

# plotting
fig = plt.figure()
ax = plt.gca()

plt.plot( xdata, derv, color = 'blue' )
plt.plot( xdata, w, color = 'red' )
plt.show()

