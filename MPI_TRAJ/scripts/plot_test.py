import matplotlib.pyplot as plt

import matplotlib.patches as mpatches
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

planck_constant = 6.62 * 10**(-34)
k = 1.23 * 10**(-23)
temperature = 300

lbs, hbs, contents = read_file( '../test', n = 3 )

means = [ 0.5 * ( lb + hb ) for lb, hb in zip( lbs, hbs )]

alpha = [ omega * j * (1 - np.exp( - omega * planck_constant / ( k * temperature ))) for omega, j in zip( means, contents) ]

lw = 2.0
plt.plot( means, alpha, linewidth = lw, color = 'k' )
plt.grid( linestyle = ':', alpha = 0.7 )

plt.show()


