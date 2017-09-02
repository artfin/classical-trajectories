import matplotlib.pyplot as plt

with open('randoms', 'r') as inputfile:
    lines = inputfile.readlines()

dots = []
for line in lines:
    if len(line) > 0:
        dots.append( float(line ))

plt.plot( dots )
plt.show()
