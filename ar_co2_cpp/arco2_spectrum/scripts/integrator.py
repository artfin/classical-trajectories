from __future__ import print_function
import numpy as np

def read_trajectory(filename):
       
    with open(filename, mode = 'r') as inputfile:
        for index, line in enumerate(inputfile):
            if index % 1000 == 0:
                print(index)

            if len(line.split()) > 0:
                if index == 0:
                    data = np.array([float(s) for s in line.split()])
                else:
                    data = np.vstack((data, np.array([float(s) for s in line.split()])))

    #R = np.array([])
    #theta = np.array([])
    #pR = np.array([])
    #pT = np.array([])
    #alpha = np.array([])
    #beta = np.array([])
    #J = np.array([])
    #phi_dot = np.array([])

    #for index, line in enumerate(lines):
        #if len(line) > 0:
           
            #if index % 1000 == 0:
                #print(index)
            
            #data = line.split()
            
            #R = np.append(R, float(data[0]))
            #theta = np.append(theta, float(data[1]))
            #pR = np.append(pR, float(data[2]))
            #pT = np.append(pT, float(data[3]))
            #alpha = np.append(alpha, float(data[4]))
            #beta = np.append(beta, float(data[5]))
            #J = np.append(J, float(data[6]))
            #phi_dot = np.append(phi_dot, float(data[7]))

    #print(R)


read_trajectory('../output/trajectory.dat')


