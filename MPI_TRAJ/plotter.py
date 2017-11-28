import sys
import matplotlib.pyplot as plt

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

    plt.xlim(( 0.0, 700.0 ))
    plt.plot( freqs, spectrum_vals, color = 'k', linewidth = 2.0, linestyle = '-' )
    plt.grid( linestyle = ':', alpha = 0.7 )
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


