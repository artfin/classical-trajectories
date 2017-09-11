import matplotlib.pyplot as plt

def read_file( filename ):
    with open( filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    x0, x1 = [], []
    for line in lines:
        
        if '#' in line:
            continue

        data = line.split()

        x0.append( float(data[0]) )
        x1.append( float(data[1]) )

    return x0, x1

x0, x1 = read_file( 'gauss_rot.txt' )

b = 200

fig, ax = plt.subplots( figsize = [8, 6] )
N, bins, patches = ax.hist( x0, bins = b, color = '#777777', normed = True )
N1, bins1, patches1 = ax.hist( x1, bins = b, color = 'red', normed = True )

plt.show()

