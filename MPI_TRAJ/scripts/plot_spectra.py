import numpy as np
import matplotlib.pyplot as plt

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
   
    freqs_range = np.zeros( len(freqs) )
    ints_range = np.zeros( len(freqs) )
    
    for i in range( int(max(freqs)) ):
        r = np.where(np.logical_and( freqs >= i, freqs <= i + 1 ))
        ints_range[i] = np.sum( np.array( ints[r] ))
        freqs_range[i] = i

    return freqs_range, ints_range
    
freqs, ints = read_file( '../output/spectra' )

freqs_range, ints_range = move_to_bins(np.array( freqs ), np.array( ints ))

plt.plot( freqs_range, ints_range, color = 'k' )


plt.show()