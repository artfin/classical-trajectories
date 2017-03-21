import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def load_data(file_path):
    with open(file_path, mode = 'r') as inputfile:
        lines = inputfile.readlines()
    
    temperatures = []
    svcs = []
    for line in lines:
        temperatures.append(int(line.split()[0]))
        svcs.append(float(line.split()[1]))
    return temperatures, svcs

def cut_data(temperatures, svcs, boundary_value):
    temperatures_temp = []
    svcs_temp = []

    for temperature, svc in zip(temperatures, svcs):
        if temperature < boundary_value:
            svcs_temp.append(svc)
            temperatures_temp.append(temperature)
        else:
            break
    return temperatures_temp, svcs_temp

parker_svc_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/parker/SVC_parker.dat'

hutson_svc_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/hutson/SVC_hutson.dat'

ab_initio_svc_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/ab-initio/SVC_ab_initio.dat'

temperatures1, parker_svc = load_data(parker_svc_path)
temperatures2, hutson_svc = load_data(hutson_svc_path)
temperatures3, ab_initio_svc = load_data(ab_initio_svc_path)

print 'len(temperatures1): {0}'.format(len(temperatures1))
print 'len(temperatures2): {0}'.format(len(temperatures2))
print 'len(temperatures3): {0}'.format(len(temperatures3))

temperatures1, parker_svc = cut_data(temperatures1, parker_svc, 300)
temperatures2, hutson_svc = cut_data(temperatures2, hutson_svc, 300)
temperatures3, ab_initio_svc = cut_data(temperatures3, ab_initio_svc, 300)

plt.plot(temperatures1, parker_svc, color = 'green')
plt.plot(temperatures2, hutson_svc, color = 'orange')
plt.plot(temperatures3, ab_initio_svc, color = 'blue')

green_patch = mpatches.Patch(color = 'green', label = 'Parker potential')
orange_patch = mpatches.Patch(color = 'orange', label = 'Huston potential')
blue_patch = mpatches.Patch(color = 'blue', label = 'New Ab-Initio potential')

plt.legend(handles = [green_patch, orange_patch, blue_patch])

plt.grid()
#plt.show()

plt.savefig('SVCs: 100-300.png')
