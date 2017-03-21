import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def load_data(filename):
    with open(filename, mode = 'r') as inputfile:
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

temperatures1, svc = load_data('SVC_hutson.dat')
temperatures2, rotational_correction = load_data('SVC1_r_hutson.dat')
temperatures3, translational_correction = load_data('SVC1_t_hutson.dat')

print 'temperatures1 len: {0}'.format(len(temperatures1))
print 'temperatures2 len: {0}'.format(len(temperatures2))
print 'temperatures3 len: {0}'.format(len(temperatures3))

boundary_temperature = 150

temperatures1, svc = cut_data(temperatures1, svc, boundary_temperature)
temperatures2, rotational_correction = cut_data(temperatures2, rotational_correction, boundary_temperature)
temperatures3, translational_correction = cut_data(temperatures3, translational_correction, boundary_temperature)

fig, ax = plt.subplots()

line1,= ax.plot(temperatures1, svc, color = 'red', label = 'SVC')
line2, = ax.plot(temperatures1, [sum(x) for x in zip(svc, rotational_correction)], color = 'blue', label = 'Rotational correction')
line3, = ax.plot(temperatures1, [sum(x) for x in zip(svc, translational_correction)], color = 'magenta', label = 'Translational correction')
line4, = ax.plot(temperatures1, [sum(x) for x in zip(svc, rotational_correction, translational_correction)], color = 'cyan', label = 'Both corrections')

line1.set_dashes([3, 1])
line2.set_dashes([3, 1])
line3.set_dashes([3, 1])
line4.set_dashes([3, 1])

red_patch = mpatches.Patch(color = 'red', label = 'SVC')
blue_patch = mpatches.Patch(color = 'blue', label = 'Rotational correction')
magenta_patch = mpatches.Patch(color = 'magenta', label = 'Translational correction')
cyan_patch = mpatches.Patch(color = 'cyan', label = 'Both corrections')

plt.legend(handles = [red_patch, blue_patch, magenta_patch, cyan_patch])

plt.grid()
plt.show()
