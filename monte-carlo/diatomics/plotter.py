import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def read_file(filename):
    with open(filename, 'r') as inputfile:
        lines = inputfile.readlines()

    temperatures = []
    constants = []
    
    for line in lines:
        temperatures.append(int(line.split()[0]))
        constants.append(float(line.split()[1]))

    return temperatures, constants

temperatures_vig, constants_vig = read_file('vigasin_diatomics.dat')
temperatures_phase, constants_phase = read_file('phase_diatomics.dat')
temperatures_rrho, constants_rrho = read_file('rrho_constant.dat')

plt.plot(temperatures_vig, constants_vig, color = 'green')
plt.plot(temperatures_phase, constants_phase, color = 'orange')
plt.plot(temperatures_rrho, constants_rrho, color = 'blue')

green_patch = mpatches.Patch(color = 'green', label = 'Vigasin formula')
orange_patch = mpatches.Patch(color = 'orange', label = 'Phase Integral')
blue_patch = mpatches.Patch(color = 'blue', label = 'RRHO')

plt.legend(handles = [green_patch, orange_patch, blue_patch])

plt.grid()
#plt.show()

plt.savefig('EquilibiumConstants.png')
