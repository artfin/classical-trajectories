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

parker_svc_path = 'parker/SVC_parker.dat'
parker_svc1_r_path = 'parker/SVC1_r_parker.dat'
parker_svc1_t_path = 'parker/SVC1_t_parker.dat'

ab_initio_svc_path = 'ab-initio/SVC_ab_initio.dat'
ab_initio_svc1_r_path = 'ab-initio/SVC1_r_ab_initio.dat'
ab_initio_svc1_t_path = 'ab-initio/SVC1_t_ab_initio.dat'

temperatures1, parker_svc = load_data(parker_svc_path)
temperatures1, parker_svc1_r = load_data(parker_svc1_r_path)
temperatures1, parker_svc1_t = load_data(parker_svc1_t_path)

print 'r: {0}'.format(len(parker_svc1_r))
print 't: {0}'.format(len(parker_svc1_t))

parker_full = [sum(x) for x in zip(parker_svc, parker_svc1_r, parker_svc1_t)]

temperatures2, ab_initio_svc = load_data(ab_initio_svc_path)
temperatures2, ab_initio_svc1_r = load_data(ab_initio_svc1_r_path)
temperatures2, ab_initio_svc1_t = load_data(ab_initio_svc1_t_path)

ab_initio_full = [sum(x) for x in zip(ab_initio_svc, ab_initio_svc1_r, ab_initio_svc1_t)]
print 'len: {0}'.format(len(ab_initio_full))
print 'len: {0}'.format(len(parker_full))
print 'len: {0}'.format(len(temperatures1))
print 'len: {0}'.format(len(temperatures2))

#temperatures1, parker_svc = cut_data(temperatures1, parker_svc, 800)
#temperatures2, hutson_svc = cut_data(temperatures2, hutson_svc, 800)
#temperatures3, ab_initio_svc = cut_data(temperatures3, ab_initio_svc, 800)

plt.plot(temperatures2, parker_full, color = 'green')
plt.plot(temperatures2, ab_initio_full, color = 'orange')

green_patch = mpatches.Patch(color = 'green', label = 'Parker w/corrections')
orange_patch = mpatches.Patch(color = 'orange', label = 'New w/corrections')

plt.legend(handles = [green_patch, orange_patch])

plt.grid()
plt.show()

#plt.savefig('SVCs: 100-800.png')
