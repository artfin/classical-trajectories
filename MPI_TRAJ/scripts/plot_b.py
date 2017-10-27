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
light_speed = 3.0 * 10**10

temperature = 300

#lbs, hbs, contents = read_file( '../spectrum300.txt', n = 3 )
#lbs, hbs, contents = read_file( '../main_test', n = 3 )
lbs, hbs, contents = read_file( '../naive_test', n = 3 )

means = [ 0.5 * ( lb + hb ) for lb, hb in zip( lbs, hbs )]

alpha = [ omega * j * (1 - np.exp( - omega * light_speed * planck_constant / ( k * temperature ))) for omega, j in zip( means, contents) ]

#N = 30 
#contents15_sm = moving_average( contents15, N )
#alpha5_sm = moving_average( alpha5, N )
#alpha10_sm = moving_average( alpha10, N )
#alpha15_sm = moving_average( alpha15, N )

lw = 2.0

begin = 10 
#plt.plot(means, contents, linewidth = lw, color = 'k' )

plt.plot( means, alpha, linewidth = lw, color = 'k' )
plt.grid( linestyle = ':', alpha = 0.7 )

plt.show()


