import matplotlib.pyplot as plt

def read_file( filename, n ):
    with open( filename, mode = 'r' ) as inputfile:
        lines = inputfile.readlines()

    lists = [ [] for _ in range(n) ]

    for line in lines:
        
        data = line.split()
        for i in range(len(data)):
            lists[i].append( float(data[i]) )

    return lists

lower_bounds, higher_bounds, contents = read_file( '../spectrum.txt', n = 3 )

means = [ ( l + h ) * 0.5 for l, h in zip( lower_bounds, higher_bounds) ]

lw = 2.0
begin = 5 
plt.plot( means[begin:], contents[begin:], linewidth = lw, color = 'k' )

plt.grid( linestyle = ':', alpha = 0.7 )
plt.show()
