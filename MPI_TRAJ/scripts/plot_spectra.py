import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

BOLTZCONST = 1.33 * 10**(-23)
PLANK_CONSTANT = 6.67 * 10**(-34)
LIGHT_SPEED = 3.00 * 10**8
Temperature = 70

CMTOHZ = 2.997 * 10**(10)

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times'

mpl.rcParams['figure.titlesize'] = 'xx-large'
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['axes.labelsize'] = 'x-large'
mpl.rcParams['axes.titlesize'] = 'large'

mpl.rcParams['xtick.labelsize'] = 'large'
mpl.rcParams['ytick.labelsize'] = 'large'

def read_file( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    freqs = []
    ints =  []
    
    for line in lines:
        data = line.split()

        freqs.append( float(data[0]) )
        ints.append( float(data[1]) )

    return freqs, ints

def move_to_bins( freqs, ints ):
   
    freqs_range = np.array([])
    ints_range = np.array([])
    
    for i in range( int(max(freqs)) ):
        r = np.where(np.logical_and( freqs >= i, freqs <= i + 1 ))
        ints_range = np.append(ints_range, np.sum( np.array( ints[r] )))
        freqs_range = np.append( freqs_range, i )

    return freqs_range, ints_range

def calc_alpha( freqs, ints ):

    alpha = np.zeros( len(freqs) )
    for i, omega, intensity in zip( range(len(freqs)), freqs, ints):
        alpha[i] = omega * (1 - np.exp( - PLANK_CONSTANT * LIGHT_SPEED * CMTOHZ * omega / (BOLTZCONST * Temperature))) * intensity
        print('omega: {0}; intensity: {1}; alpha: {2}'.format(omega, intensity, alpha[i]))
    
    return alpha


freqs, ints = read_file( '../first_exp/spectra_weighted' )

freqs_range, ints_range = move_to_bins(np.array( freqs ), np.array( ints ))

# ----------------------------------------------
# absorption coefficient
#alpha = calc_alpha( freqs_range, ints_range )

#plt.plot( freqs_range, alpha, color = 'k' )
#plt.show()
# ----------------------------------------------


# ----------------------------------------------
# Spectral function
plt.plot( freqs_range[10:], ints_range[10:], color = 'k' )

#plt.title(r"\Large Spectral function $J(\omega)$")
#plt.semilogy( freqs_range, ints_range, color = 'k' )

#plt.xlim( (0, 170) )
#plt.ylim( (0, 0.005) )

#plt.xlabel(r"$\omega$, \textbf{cm}$^{-1}$")

plt.grid( linestyle = ':', alpha = 0.7 )
plt.show()
# ----------------------------------------------


