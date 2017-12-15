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

def plot_spectrum( filename ):
    
    freqs, spectrum_vals = read_file( filename )
    spectrum_vals = [ v * 10**7 for v in spectrum_vals ]
    spectrum_vals_corr = [ v * 0.963 for v in spectrum_vals ]

    freqs_bg, spectrum_bg = read_file( './pictures/bosomworth_gash1965_295K.txt' )
    freqs_kw, spectrum_kw = read_file( './pictures/kiss_welsh.txt' )
    freqs_mq, spectrum_mq = read_file( './pictures/mcquarrie1968_295K.txt' )

    freqs_matlab_buryak_fernandez, spectrum_matlab_buryak_fernandez = read_file( './pictures/buryak_matlab_spectrum_fernandez_295K.txt' )
    spectrum_matlab_buryak_fernandez = norm( spectrum_matlab_buryak_fernandez )
    
    freqs_matlab_buryak_meyer, spectrum_matlab_buryak_meyer = read_file( './pictures/buryak_matlab_spectrum_meyer_295K.txt' )
    spectrum_matlab_buryak_meyer = norm( spectrum_matlab_buryak_meyer )

    fig = plt.figure()
    plt.title( r'\LARGE He-Ar binary adsorption coefficient, 295 K')
    plt.xlabel(r'$\omega$, cm$^{-1}$') 
    plt.ylabel(r'$\alpha$, cm$^{-1} \cdot$ amagat$^{-2}$')

    plt.xlim(( 0.0, 700.0 ))
    l1, = plt.plot( freqs, spectrum_vals, color = 'k', linewidth = 2.0, linestyle = ':' )
    l2 = plt.scatter( freqs_bg, spectrum_bg, color = 'k', s = 20 )
    l3, = plt.plot( freqs, spectrum_vals_corr, color = 'k', lw = 2.0, linestyle = '-' )
    l5, = plt.plot( freqs_matlab_buryak_fernandez, spectrum_matlab_buryak_fernandez, color = 'r' )
    l6, = plt.plot( freqs_matlab_buryak_meyer, spectrum_matlab_buryak_meyer, color = 'b' )
    plt.grid( linestyle = ':', alpha = 0.7 )
    
    fig.legend((l1, l2, l3, l5, l6), ( r'My MPI program', r'Bosomworth, Gass, 1965', r'My data multiplied', r'I. Buryak (Fernandez)', r'I. Buryak (Meyer)'), 'lower center', ncol = 4, fancybox = True, shadow = True)

    plt.tight_layout()
    plt.show()


if ( len(sys.argv) != 3 ):
    print("Need to call with 2 arguments: ")
    print("1: path-to-file")
    print("2: type of plot: specfunc or spectrum")
    sys.exit( 1 )

if ( sys.argv[2] == "specfunc" ):
    plot_specfunc( sys.argv[1] )

if ( sys.argv[2] == "spectrum" ):
    plot_spectrum( sys.argv[1] )


