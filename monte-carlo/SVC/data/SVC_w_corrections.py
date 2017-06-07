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

#parker_svc_path = 'parker/SVC_parker.dat'
#parker_svc1_r_path = 'parker/SVC1_r_parker.dat'
#parker_svc1_t_path = 'parker/SVC1_t_parker.dat'

ab_initio_svc_path = 'ab-initio/SVC_from0.dat'
ab_initio_svc1_r_path = 'ab-initio/SVC1r_long.dat'
ab_initio_svc1_t_path = 'ab-initio/SVC1t_long.dat'

#temperatures1, parker_svc = load_data(parker_svc_path)
#temperatures1, parker_svc1_r = load_data(parker_svc1_r_path)
#temperatures1, parker_svc1_t = load_data(parker_svc1_t_path)

#print 'r: {0}'.format(len(parker_svc1_r))
#print 't: {0}'.format(len(parker_svc1_t))

#parker_full = [sum(x) for x in zip(parker_svc, parker_svc1_r, parker_svc1_t)]

temperatures1, ab_initio_svc = load_data(ab_initio_svc_path)
temperatures2, ab_initio_svc1_r = load_data(ab_initio_svc1_r_path)
temperatures3, ab_initio_svc1_t = load_data(ab_initio_svc1_t_path)

ub = 200

temperatures1 = temperatures1[:ub]
ab_initio_svc = ab_initio_svc[:ub]

temperatures2 = temperatures2[:ub]
ab_initio_svc1_r = ab_initio_svc1_r[:ub]

temperatures_temp = temperatures3[:ub]
ab_initio_svc1_t = ab_initio_svc1_t[:ub]

ab_initio_full = [sum(x) for x in zip(ab_initio_svc, ab_initio_svc1_r, ab_initio_svc1_t)]

res = [y-x for x, y in zip(ab_initio_svc, ab_initio_full)]

print(ab_initio_full[0] - ab_initio_svc[0])

fig = plt.figure()
plt.rc('text', usetex = True)
lw = 1.75

l, = plt.plot(temperatures1, res, color = '0.5', linestyle = 'solid', linewidth = lw)
#l1, = plt.plot(temperatures1, ab_initio_svc, color = '0.5', linestyle = 'dashed', linewidth = lw) 
#l2, = plt.plot(temperatures1, ab_initio_full, color = '0.6', linestyle = 'solid', linewidth = lw)

plt.xlabel(r'\textbf{T}, (K)')
plt.ylabel(r'Quantum corrections to B$_2$, $\left( \, {\Large\displaystyle\frac{\textit{cm}^{\, 3}}{\textit{mol}}} \, \right)$')

#fig.legend((l1, l2), (r'SVC, \textit{ab initio}', r'SVC w/corrections, \textit{ab initio}'), 'lower center', ncol = 2, fancybox = True, shadow = True, prop = {'size': 'large'})
fig.legend((l, ), (r'SVC quantum corrections, \textit{ab initio} potential', ), 'lower center', ncol = 1, fancybox = True, shadow = True, prop = {'size': 'large'})

plt.grid(linestyle = ':', alpha = 0.7)
plt.show()

#plt.savefig('SVCs: 100-800.png')
