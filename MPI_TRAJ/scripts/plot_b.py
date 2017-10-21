import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.interpolate import spline
import numpy as np

def read_file( filename, n ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    lists = [ [] for _ in range(n) ]

    for line in lines:
        
        data = line.split()
        for i in range(len(data)):
            lists[i].append( float(data[i]) )

    return lists

def moving_average( arr, N ):

    cumsum, moving_aves = [0], []

    for i, x in enumerate(arr, 1):
        cumsum.append( cumsum[i-1] + x )
        if i >= N:
            moving_ave = ( cumsum[i] - cumsum[i - N]) / N
            moving_aves.append( moving_ave )
    
    return moving_aves

planck_constant = 6.62 * 10**(-34)
k = 1.23 * 10**(-23)
temperature = 300

#lbs1, hbs1, contents1 = read_file( '../experiments/b_selection/spectrum_b_1.txt', n = 3 )
#lbs5, hbs5, contents5 = read_file( '../experiments/b_selection/spectrum_b_5.txt', n = 3 )
#lbs10, hbs10, contents10 = read_file( '../experiments/b_selection/spectrum_b_10.txt', n = 3 )
lbs15, hbs15, contents15 = read_file( '../experiments/b_selection/spectrum_b_15.txt', n = 3 )

#means1 = [ 306.0 / 600.0 * 0.5 * ( lb + hb ) for lb, hb in zip( lbs1, hbs1 )]
#means5 = [ 306.0 / 600.0 * 0.5 * ( lb + hb ) for lb, hb in zip( lbs5, hbs5 )]
#means10 = [ 306.0 / 600.0 * 0.5 * ( lb + hb ) for lb, hb in zip( lbs10, hbs10 )]
means15 = [ 0.5 * ( lb + hb ) for lb, hb in zip( lbs15, hbs15 )]

#alpha1 = [ omega * j * (1 - np.exp( - omega * planck_constant / ( k * temperature ))) for omega, j in zip( means1, contents1) ]
#alpha5 = [ omega * j * (1 - np.exp( - omega * planck_constant / ( k * temperature ))) for omega, j in zip( means5, contents5) ]
#alpha10 = [ omega * j * (1 - np.exp( - omega * planck_constant / ( k * temperature ))) for omega, j in zip( means10, contents10) ]
alpha15 = [ omega * j * (1 - np.exp( - omega * planck_constant / ( k * temperature ))) for omega, j in zip( means15, contents15) ]

#N = 30 
#contents15_sm = moving_average( contents15, N )
#alpha5_sm = moving_average( alpha5, N )
#alpha10_sm = moving_average( alpha10, N )
#alpha15_sm = moving_average( alpha15, N )

lw = 2.0

begin = 10 
#plt.plot(means15[begin:(len(means15) - N + 1)], contents15_sm[begin:], linewidth = lw, color = 'k' )

#plt.plot( means5[:(len(means5) - N + 1)], alpha5_sm, linewidth = lw, color = 'g' )
#plt.plot( means10[:(len(means10) - N + 1)], alpha10_sm, linewidth = lw, color = 'r' )
#plt.plot( means15[:(len(means15) - N + 1)], alpha15_sm, linewidth = lw, color = 'k' )
plt.plot( means15[begin:], contents15[begin:], linewidth = lw, color = 'k' )
plt.grid( linestyle = ':', alpha = 0.7 )

plt.show()


