import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math as m

def load_data(filename):
	with open(filename, mode = 'r') as inputfile:
		data = inputfile.readlines()

	temperatures = []
	constants = []

	for line in data:
		_ = line.split()
		temperatures.append(1. / float(_[0]))
		constants.append(m.log(float(_[1])))

	return temperatures, constants

temperatures, parker = load_data('parker.dat')
temperatures, hutson = load_data('hutson.dat')

temperatures = temperatures[10:]
parker = parker[10:]
hutson = hutson[10:]

plt.plot(temperatures, parker, 'red')
plt.plot(temperatures, hutson, 'green')

patch1 = mpatches.Patch(color = 'red', label = 'Parker-Snow')
patch2 = mpatches.Patch(color = 'green', label = 'Hutson')
plt.legend(handles = [patch1, patch2])
plt.grid()
# plt.savefig('Both.png')
plt.show()
