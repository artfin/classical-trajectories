import sys
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.unicode'] = True

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times'

mpl.rcParams['figure.titlesize'] = "xx-large"
mpl.rcParams['legend.fontsize'] = "large"
mpl.rcParams['axes.labelsize'] = "x-large"
mpl.rcParams['axes.titlesize'] = "large"

mpl.rcParams['xtick.labelsize'] = "large"
mpl.rcParams['ytick.labelsize'] = "large"

def norm( arr ):
    return [ v * 10**7 for v in arr ]

def read_file( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    freqs, intensities = [], []

    for line in lines:
        data = line.split()

        freqs.append( float(data[0]) )
        intensities.append( float(data[1]) )

    return freqs, intensities

def plot_specfunc( filename ):
    freqs, specfunc_vals = read_file( filename )

    maximum_value = max( specfunc_vals )

    plt.gca().set_yscale('log')
    plt.xlim(( 0.0, 500.0 ))
    plt.ylim(( maximum_value * 1.10 * 0.01, maximum_value * 1.10))

    plt.plot( freqs, specfunc_vals, color = 'k', linewidth = 2.0, linestyle = '-' )
    plt.grid( linestyle = ':', alpha = 0.7 )
    plt.show()

def plot_spectrum():
    
    freqs_class, spectrum_class = read_file( './spectrum_buryak_dipole.txt' )
    spectrum_class = norm( spectrum_class )

    freqs_d1, spectrum_d1 = read_file( './spectrum_buryak_dipole_d1.txt' )
    spectrum_d1 = norm( spectrum_d1 )

    freqs_d2, spectrum_d2 = read_file( './spectrum_buryak_dipole_d2.txt' )
    spectrum_d2 = norm( spectrum_d2 )

    freqs_d3, spectrum_d3 = read_file( './spectrum_buryak_dipole_d3.txt' )
    spectrum_d3 = norm( spectrum_d3 )

    freqs_bg, spectrum_bg = read_file( '../pictures/bosomworth_gash1965_295K.txt' )

    fig = plt.figure()
    plt.title( r'\LARGE He-Ar binary adsorption coefficient, 295 K')
    plt.xlabel(r'$\omega$, cm$^{-1}$') 
    plt.ylabel(r'$\alpha \cdot 10^{7}$, cm$^{-1} \cdot$ amagat$^{-2}$')

    plt.xlim(( 0.0, 700.0 ))
    l1, = plt.plot( freqs_class, spectrum_class, color = 'k', linewidth = 2.0, linestyle = '-' )
    l2 = plt.scatter( freqs_bg, spectrum_bg, color = 'k', s = 20 )
    l3, = plt.plot( freqs_d1, spectrum_d1, color = 'r', linewidth = 2.0, linestyle = '-') 
    l4, = plt.plot( freqs_d2, spectrum_d2, color = 'b', linewidth = 2.0, linestyle = '-')
    l5, = plt.plot( freqs_d3, spectrum_d3, color = 'g', linewidth = 2.0, linestyle = '-')
    plt.grid( linestyle = ':', alpha = 0.7 )
    
    fig.legend((l1, l2, l3, l4, l5), ( r'Class.', r'Bosomworth, Gass, 1965', 'D1', 'D2', 'D3'), 'lower center', ncol = 5, fancybox = True, shadow = True)

    plt.tight_layout()
    plt.show()


if ( len(sys.argv) != 2 ):
    print("1: type of plot: specfunc or spectrum")
    sys.exit( 1 )

if ( sys.argv[1] == "specfunc" ):
    plot_specfunc()

if ( sys.argv[1] == "spectrum" ):
    plot_spectrum()


