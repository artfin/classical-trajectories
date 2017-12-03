def read_file( filename ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    words = lines[0].split()

    arr = []
    for word in words:
        arr.append( float(word) )

    return arr

freqs = read_file( './freqs.txt' )
alphas = read_file( './spectrum.txt' )

with open( '1.txt', mode = 'w' ) as out:
    for freq, alpha in zip( freqs, alphas ):
        out.write( str(freq) + ' ' + str(alpha) + '\n' )
