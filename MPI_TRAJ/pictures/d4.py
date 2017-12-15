import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as inter

from pprint import pprint

def read_file( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    freqs, specfunc = [], []
    for line in lines:
        if len(line) > 0:
            data = line.split()

        freqs.append( float(data[0]) )
        specfunc.append( float(data[1]) )

    return freqs, specfunc

def symmetrize( freqs, specfunc ):
    freqs = freqs + [ -f for f in reversed(freqs) ][:-1]
    specfunc = specfunc + [ s for s in reversed(specfunc) ][:-1]
    return freqs, specfunc

freqs, specfunc = read_file( "./specfunc_buryak_fit_dipole.txt" )

freqs, specfunc = symmetrize( freqs, specfunc )

freqs_delta = (freqs[1] - freqs[0]) * 3.0 * 10**10
print('freqs_delta: ' + str(freqs_delta))

time_delta = 1.0 / ( len(freqs) * freqs_delta )
print('time_delta: ' + str(time_delta))

time = np.linspace( -len(freqs) / 2.0 * time_delta,
                     len(freqs) / 2.0 * time_delta,
                     len(freqs) )

time_shifted = np.fft.ifftshift( time )
specfunc_ifft_shifted = np.fft.ifft( specfunc )
specfunc_ifft = np.fft.ifftshift( specfunc_ifft_shifted )

#plt.scatter( freqs, specfunc, s = 2 )
#plt.ylim( (0, 6e-76) )
#plt.show()
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.set_yscale('log')

pairs = []
for el1, el2 in zip(time, specfunc_ifft.real ):
    if el1 > 0:
        pairs.append( (el1, el2) )

spl = inter.InterpolatedUnivariateSpline( 
        [ p[0] for p in pairs], 
        [ p[1].real for p in pairs ]
    )

for t, v in zip(time, specfunc_ifft):
    print('t: {0}; spl(t): {1}; v: {2}'.format(t, spl(t), v))

x = np.linspace( 0, time[-1], 500 ) 

plt.xlim( (0, 2e-12) )
plt.scatter( time, specfunc_ifft.real, s = 3 )
plt.plot( x, spl(x), color = 'k', linestyle = ':', lw = 2 )
plt.ylim( (-5.0e-80, 1.0e-79) )
plt.show()

