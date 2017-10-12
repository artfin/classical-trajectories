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
temperature = 100

lower_bounds_20, higher_bounds_20, contents_20 = read_file( '../experiments/spectrum_lconst_20.txt', n = 3 )
lower_bounds_50, higher_bounds_50, contents_50 = read_file( '../experiments/spectrum_lconst_50.txt', n = 3 )
lower_bounds_100, higher_bounds_100, contents_100 = read_file( '../experiments/spectrum_lconst_100.txt', n = 3 )

means_20 = [ ( l + h ) * 0.5 for l, h in zip( lower_bounds_20, higher_bounds_20) ]
means_50 = [ ( l + h ) * 0.5 for l, h in zip( lower_bounds_50, higher_bounds_50) ]
means_100 = [ ( l + h ) * 0.5 for l, h in zip( lower_bounds_100, higher_bounds_100) ]

alpha_20 = [ omega**2 * j * (1 - np.exp( - omega * planck_constant / ( k * temperature ))) for omega, j in zip( means_20, contents_20) ]
alpha_50 = [ omega**2 * j * (1 - np.exp( - omega * planck_constant / ( k * temperature ))) for omega, j in zip( means_50, contents_50) ]
alpha_100 = [ omega**2 * j * (1 - np.exp( - omega * planck_constant / ( k * temperature ))) for omega, j in zip( means_100, contents_100) ]

lw = 2.0
begin = 0
end = 20 
plt.plot( means_20[begin:end], contents_20[begin:end], linewidth = lw, color = 'r' )
plt.plot( means_50[begin:end], contents_50[begin:end], linewidth = lw, color = 'g' )
plt.plot( means_100[begin:end], contents_100[begin:end], linewidth = lw, color = 'b' )

patch_red = mpatches.Patch( color = 'r', label = r'$L^2$ = 20' )
patch_green = mpatches.Patch( color = 'g', label = r'$L^2$ = 50' )
patch_blue = mpatches.Patch( color = 'b', label = r'$L^2$ = 100' )

#plt.plot( means_20[begin:end], alpha_20[begin:end], linewidth = lw, color = 'r' )
#plt.plot( means_50[begin:end], alpha_50[begin:end], linewidth = lw, color = 'g' )
#plt.plot( means_100[begin:end], alpha_100[begin:end], linewidth = lw, color = 'b' )
plt.grid( linestyle = ':', alpha = 0.7 )
plt.legend( handles = [patch_red, patch_green, patch_blue] )
plt.show()
