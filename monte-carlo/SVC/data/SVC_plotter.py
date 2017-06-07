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
ab_initio_svc1_r_path = 'ab-initio/SVC1r_long.dat'
ab_initio_svc1_t_path = 'ab-initio/SVC1t_long.dat'

temperatures1, parker_svc = load_data(parker_svc_path)
temperatures2, hutson_svc = load_data(hutson_svc_path)
temperatures3, ab_initio_svc = load_data(ab_initio_svc_path)
temperatures4, ab_initio_svc1_r = load_data(ab_initio_svc1_r_path)
temperatures5, ab_initio_svc1_t = load_data(ab_initio_svc1_t_path)

lowest_temperature = 200
highest_temperature = 500
temperatures1, parker_svc = cut_data(temperatures1, parker_svc, lowest_temperature, highest_temperature)
temperatures2, hutson_svc = cut_data(temperatures2, hutson_svc, lowest_temperature, highest_temperature)
temperatures3, ab_initio_svc = cut_data(temperatures3, ab_initio_svc, lowest_temperature, highest_temperature)
temperatures4, ab_initio_svc1_r = cut_data(temperatures4, ab_initio_svc1_r, lowest_temperature, highest_temperature)
temperatures5, ab_initio_svc1_t = cut_data(temperatures5, ab_initio_svc1_t, lowest_temperature, highest_temperature)

print(ab_initio_svc1_r)

def plot_svcs():
    fig = plt.figure()

    plt.rc('text', usetex = True)
    lw = 1.75
    
    l1, = plt.plot(temperatures1, parker_svc, color = '0.6', linestyle = 'solid', linewidth = lw)
    l2, = plt.plot(temperatures2, ab_initio_svc, color = '0.5', linestyle = 'dashed', linewidth = lw)
    l3, = plt.plot(temperatures3, hutson_svc, color = '0.4', linestyle = 'dotted', linewidth = lw)

    plt.xlabel(r'\textbf{T}, (K)')
    plt.ylabel(r'$\textbf{B}_2$, $\left( \, {\Large\displaystyle\frac{\textit{cm}^{\, 3}}{\textit{mol}}} \, \right)$')

    #fig.legend((l1, l2, l3), ('Parker', 'Ab-Initio', 'Hutson'), 'lower center', ncol = 3, fancybox = True, shadow = True, prop = {'size': 'large'})
    #plt.grid()
    #plt.show()

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

    for temperature, svc, error in zip(temperatures, svcs, errors):
        l4 = plt.scatter(temperature, svc, marker = '*', color = 'k') 
        plt.errorbar(temperature, svc, error, color = '0.3', fmt = "*", capsize = 3)

    fig.legend((l1, l2, l3, l4), ('PSP', r'\textit{ab initio}', 'Hutson', 'Experiment'), 'lower center', ncol = 4, fancybox = True, shadow = True, prop = {'size': 'large'})

    plt.grid(linestyle = ':', alpha = 0.7)
    plt.show()

temperatures = [213., 213., 223., 223.2, 242., 242., 248.2, 262., 273.2, 276., 276., 288.2, 
        290.0, 290.0, 295., 295., 296., 296.15, 300.,
        300.0, 303.15, 303.2, 310., 313.2, 320., 322.85, 323.1, 330., 330., 333.15, 363.15,
        365., 400., 400., 425., 450., 450., 475.]

svcs = [-86.3, -94.0, -75.5, -74.8, -62.9, -70.0, -58.4, -50.8, -50.6, -43.4, -51.0, -40.3,
        -45.2, -46.4, -37.2, -44.0, -37.0, -44.1, -40.8,
        -41.7, -31.8, -34.2, -38.6, -31.2, -35.3, -30.1, -28.3, -27.3, -35.0, -25.8, -19.6,
        -16.2, -6.0, -13.0, -3.1, 0.5, -7.0, 1.7]

errors = [5., 7., 5.0, 1.0, 5.0, 7.0, 1.0, 5.0, 1.0, 5.0, 6.0, 2.0,
          1.4, 4.0, 5.0, 6.0, 2.0, 5.0, 1.3,
          4.0, 4.6, 2.0, 4.0, 2.0, 1.3, 2.0, 2.0, 5.0, 5.0, 4.2, 4.2,
          5.0, 4.0, 3.0, 4.0, 4.0, 2.0, 4.0]


plot_svcs()

def plot_diff_ab():
    fig = plt.figure()

    plt.rc('text', usetex = True)
    lw = 1.75
    
    plt.plot((200.0, 500.0), (0.0, 0.0), color = 'k', linestyle = 'dotted', linewidth = 2)

    summab, summhut = 0, 0
    count = 0

    for _t, _s, _e in zip(temperatures, svcs, errors):
        
        t = int(_t)
        index = t - 200 - 1
    
        dab = ab_initio_svc[index] - _s
        dhut = hutson_svc[index] - _s

        #print('index: {0}; t: {1}; svc exp: {2}; svc ab initio: {3}; d: {4}'.format(index, t, _s, ab_initio_svc[index], d))
        
        l = plt.scatter(_t, dab, color = 'k', marker = '*')
        plt.errorbar(_t, dab, _e, color = 'k', capsize = 2)
        
        summab += dab**2
        summhut += dhut**2

        count += 1

    print('ab initio: {0}; hutson: {1}'.format(summab / count, summhut / count))
    
    plt.xlabel(r'\textbf{T}, (K)')
    plt.ylabel(r'$\Delta \textbf{B}_2$, $\left( \, {\Large\displaystyle\frac{\textit{cm}^{\, 3}}{\textit{mol}}} \, \right)$')
    
    fig.legend((l,), (r'$\Delta$ (SVC (\textit{ab initio}), SVC (\textit{experiment}))', ), 'lower center', ncol = 1, fancybox = True, shadow = True, prop = {'size': 'large'})

    plt.text(475.0, 3.0, r'${\Huge \sigma^2 = 16.319}$', color = 'k', bbox = dict(facecolor='none', edgecolor = 'black', boxstyle = 'round,pad=1'))

    plt.grid(linestyle = ':', alpha = 0.7)
    plt.show()

def plot_diff_hut():
    fig = plt.figure()

    plt.rc('text', usetex = True)
    lw = 1.75
    
    plt.plot((200.0, 500.0), (0.0, 0.0), color = 'k', linestyle = 'dotted', linewidth = 2)

    summab, summhut = 0, 0
    count = 0

    for _t, _s, _e in zip(temperatures, svcs, errors):
        
        t = int(_t)
        index = t - 200 - 1
    
        dab = ab_initio_svc[index] - _s
        dhut = hutson_svc[index] - _s

        #print('index: {0}; t: {1}; svc exp: {2}; svc ab initio: {3}; d: {4}'.format(index, t, _s, ab_initio_svc[index], d))
        
        l = plt.scatter(_t, dhut, color = 'k', marker = '*')
        plt.errorbar(_t, dhut, _e, color = 'k', capsize = 2)
        
        summab += dab**2
        summhut += dhut**2

        count += 1

    print('ab initio: {0}; hutson: {1}'.format(summab / count, summhut / count))
    
    plt.xlabel(r'\textbf{T}, (K)')
    plt.ylabel(r'$\Delta \textbf{B}_2$, $\left( \, {\Large\displaystyle\frac{\textit{cm}^{\, 3}}{\textit{mol}}} \, \right)$')
    
    fig.legend((l,), (r'$\Delta$ (SVC (\textit{hutson}), SVC (\textit{experiment}))', ), 'lower center', ncol = 1, fancybox = True, shadow = True, prop = {'size': 'large'})

    plt.text(475.0, 3.0, r'${\Huge \sigma^2 = 21.552}$', color = 'k', bbox = dict(facecolor='none', edgecolor = 'black', boxstyle = 'round,pad=1'))

    plt.grid(linestyle = ':', alpha = 0.7)
    plt.show()

def plot_diff_ab_wcor():
    fig = plt.figure()

    plt.rc('text', usetex = True)
    lw = 1.75
    
    plt.plot((200.0, 500.0), (0.0, 0.0), color = 'k', linestyle = 'dotted', linewidth = 2)

    summab = 0
    count = 0
    
    for _t, _s, _e in zip(temperatures, svcs, errors):
        
        t = int(_t)
        index = t - 200 - 1
        
        if index < len(ab_initio_svc1_r):
            dab = ab_initio_svc[index] + ab_initio_svc1_r[index] + ab_initio_svc1_t[index] - _s
            print('index: {0}; svc: {1}; svc1r: {2}; svc1t: {3}'.format(index, ab_initio_svc[index], ab_initio_svc1_r[index], ab_initio_svc1_t[index]))
        else:
            ab = ab_initio_svc[index] - _s

        l = plt.scatter(_t, dab, color = 'k', marker = '*')
        plt.errorbar(_t, dab, _e, color = 'k', capsize = 2)
        
        summab += dab**2

        count += 1

    print('ab initio: {0}'.format(summab / count))
    
    plt.xlabel(r'\textbf{T}, (K)')
    plt.ylabel(r'$\Delta \textbf{B}_2$, $\left( \, {\Large\displaystyle\frac{\textit{cm}^{\, 3}}{\textit{mol}}} \, \right)$')
    
    fig.legend((l,), (r'$\Delta$ (SVC (\textit{ab initio}), SVC (\textit{experiment}))', ), 'lower center', ncol = 1, fancybox = True, shadow = True, prop = {'size': 'large'})

    plt.text(475.0, 5.0, r'${\Huge \sigma^2 = 13.487}$', color = 'k', bbox = dict(facecolor='none', edgecolor = 'black', boxstyle = 'round,pad=1'))

    plt.grid(linestyle = ':', alpha = 0.7)
    plt.show()

#plot_diff_ab_wcor()
