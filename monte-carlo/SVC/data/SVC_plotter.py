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

parker_svc_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/parker/SVC_parker.dat'
hutson_svc_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/hutson/SVC_from0.dat'
ab_initio_svc_path = '/home/artfin/Desktop/repos/classical-trajectories/classical-trajectories/monte-carlo/SVC/data/ab-initio/SVC_from0.dat'

temperatures1, parker_svc = load_data(parker_svc_path)
temperatures2, hutson_svc = load_data(hutson_svc_path)
temperatures3, ab_initio_svc = load_data(ab_initio_svc_path)

print 'len(temperatures1): {0}'.format(len(temperatures1))
print 'len(temperatures2): {0}'.format(len(temperatures2))
print 'len(temperatures3): {0}'.format(len(temperatures3))


lowest_temperature = 200
highest_temperature = 500

temperatures1, parker_svc = cut_data(temperatures1, parker_svc, lowest_temperature, highest_temperature)
temperatures2, hutson_svc = cut_data(temperatures2, hutson_svc, lowest_temperature, highest_temperature)
temperatures3, ab_initio_svc = cut_data(temperatures3, ab_initio_svc, lowest_temperature, highest_temperature)

plt.plot(temperatures1, parker_svc, color = 'green', linestyle = 'dashed')
plt.plot(temperatures2, hutson_svc, color = 'red', linestyle = 'dashdot')
plt.plot(temperatures3, ab_initio_svc, color = 'blue', linestyle = 'dotted')

green_patch = mpatches.Patch(color = 'green', label = 'Parker potential')
orange_patch = mpatches.Patch(color = 'red', label = 'Hutson potential')
blue_patch = mpatches.Patch(color = 'blue', label = 'Ab-Initio potential')


#temperatures = [213., 223., 242., 262., 276., 288.2, 296., 303.2, 313.2, 323.1, 333.2, 365.]
#single_rep = [-93.1, -84.8, -71.3, -59.8, -53.0, -47.8, -44.7, -42.0, -38.5, -35.4, -32.4, -24.2]
#split_rep = [-88.3, -80.2, -67.3, -56.3, -49.7, -44.7, -41.7, -39.1, -35.8, -32.8, -29.9, -22.0] 
#experiment = [-93.9, -84.5, -67.2, -53.8, -48.7, -40.3, -37.05, -34.21, -31.20, -28.30, -25.8, -19.6] 
#errors = [10., 10., 10., 10., 7., 7., 7., 7., 7., 7., 7., 7.]

#for temperature, rep1, rep2, svc, error in zip(temperatures, single_rep, split_rep, experiment, errors):
    #plt.scatter(temperature, svc, marker = '*', color = 'k')
    #plt.scatter(temperature, rep1, marker = 'x', color = 'k')
    #plt.scatter(temperature, rep2, marker = 'x', color = 'magenta')
    #plt.errorbar(temperature, svc, error, color = 'red')

#temperatures = [213., 213., 223., 223.2, 242., 242., 248.2, 262., 273.2, 276., 276., 288.2, 
        #290.0, 290.0, 295., 295., 296., 296.15, 300.,
        #300.0, 303.15, 303.2, 310., 313.2, 320., 322.85, 323.1, 330., 330., 333.15, 363.15,
        #365., 400., 400., 425., 450., 450., 475.]

#svcs = [-86.3, -94.0, -75.5, -74.8, -62.9, -70.0, -58.4, -50.8, -50.6, -43.4, -51.0, -40.3,
        #-45.2, -46.4, -37.2, -44.0, -37.0, -44.1, -40.8,
        #-41.7, -31.8, -34.2, -38.6, -31.2, -35.3, -30.1, -28.3, -27.3, -35.0, -25.8, -19.6,
        #-16.2, -6.0, -13.0, -3.1, 0.5, -7.0, 1.7]

#errors = [5., 7., 5.0, 1.0, 5.0, 7.0, 1.0, 5.0, 1.0, 5.0, 6.0, 2.0,
          #1.4, 4.0, 5.0, 6.0, 2.0, 5.0, 1.3,
          #4.0, 4.6, 2.0, 4.0, 2.0, 1.3, 2.0, 2.0, 5.0, 5.0, 4.2, 4.2,
          #5.0, 4.0, 3.0, 4.0, 4.0, 2.0, 4.0]


#for temperature, svc, error in zip(temperatures, svcs, errors):
    #plt.scatter(temperature, svc, marker = '*', color = 'k') 
    ##plt.errorbar(temperature, svc, error, color = 'k', fmt = "*--")

#error_patch = mpatches.Patch(color = 'k', label = 'Dymond-Smith Data')

#plt.legend(handles = [error_patch, orange_patch, blue_patch])
plt.legend(handles = [green_patch, orange_patch, blue_patch])

plt.grid()
plt.show()

#plt.savefig('Exp_Hutson_Julia1.png')
