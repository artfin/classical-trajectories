import matplotlib.pyplot as plt
import math as m

def read_file( filename ):
    with open(filename, 'r') as inputfile:
        lines = inputfile.readlines()

    jx = []
    jy = []
    pR = []

    for line in lines:
        data = line.split()
        jx.append( float(data[0]) )
        jy.append( float(data[1]) )
        pR.append( float(data[2]) )

    return jx, jy, pR

jx1, jy1, pR1 = read_file('ic.txt')
jx2, jy2, pR2 = read_file('ic2.txt')

fig, ax = plt.subplots(figsize=[8,6])

N1, bins1, patches1 = ax.hist(pR1, bins = 100, normed = True )
N2, bins2, patches2 = ax.hist(pR2, bins = 100, normed = True )

# plt.plot(jx, color = 'k', linewidth = 1.2)
# plt.plot(jx, color = 'k', linewidth = 0.5 )
plt.show()
