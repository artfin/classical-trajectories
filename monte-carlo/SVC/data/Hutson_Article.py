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

def cut_data(temperatures, svcs, low_bound, high_bound):
    temperatures_temp = []
    svcs_temp = []

    for temperature, svc in zip(temperatures, svcs):
        if temperature < high_bound and temperature > low_bound:
            svcs_temp.append(svc)
            temperatures_temp.append(temperature)
    return temperatures_temp, svcs_temp


hutson_svc_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/hutson/SVC_from0.dat'

temperatures, hutson_svc = load_data(hutson_svc_path)

lowest_temperature = 200
highest_temperature = 500

temperatures, hutson_svc = cut_data(temperatures, hutson_svc, lowest_temperature, highest_temperature)

plt.plot(temperatures, hutson_svc, color = 'orange')

orange_patch = mpatches.Patch(color = 'orange', label = 'Hutson potential')


temperatures = [213., 223., 242., 262., 276., 288.2, 296., 303.2, 313.2, 323.1, 333.2, 365.]
single_rep = [-93.1, -84.8, -71.3, -59.8, -53.0, -47.8, -44.7, -42.0, -38.5, -35.4, -32.4, -24.2]
split_rep = [-88.3, -80.2, -67.3, -56.3, -49.7, -44.7, -41.7, -39.1, -35.8, -32.8, -29.9, -22.0] 
experiment = [-93.9, -84.5, -67.2, -53.8, -48.7, -40.3, -37.05, -34.21, -31.20, -28.30, -25.8, -19.6] 
errors = [10., 10., 10., 10., 7., 7., 7., 7., 7., 7., 7., 7.]

for temperature, rep1, rep2, svc, error in zip(temperatures, single_rep, split_rep, experiment, errors):
    plt.scatter(temperature, svc, marker = '*', color = 'k')
    plt.scatter(temperature, rep1, marker = 'x', color = '#4f6d27')
    plt.scatter(temperature, rep2, marker = 'x', color = 'red')
    plt.errorbar(temperature, svc, error, color = 'k', fmt = '*--')

error_patch = mpatches.Patch(color = 'k', label = 'Hutson Experimental Data')

green_patch = mpatches.Patch(color = '#4f6d27', label = 'Single Repulsion Potential')
red_patch = mpatches.Patch(color = 'red', label = 'Split Repulsion Potential')

plt.legend(handles = [orange_patch, error_patch, green_patch, red_patch])

plt.grid()
plt.show()

#plt.savefig('Hutson_article.png')
