import numpy as np
import matplotlib.pyplot as plt

temperatures = [213., 213., 223., 233., 242., 242., 248.2, 262., 273.2, 276., 276., 288.2,
        290.0, 290.0, 295.0, 295.0, 296.0, 296.15, 300.,
        300.0, 303.15, 303.2, 310.0, 313.2, 320.0, 322.85, 323.1, 330.0, 330.0, 333.15, 363.15,
        365.0, 400.0, 400.0, 425., 450.0, 450.0, 475.]
svcs = [-86.3, -94.0, -75.5, -74.8, -62.9, -70.0, -58.4, -50.8, -50.6, -43.4, -51.0, -40.3,
        -45.2, -46.4, -37.2, -44.0, -37.0, -44.1, -40.8,
        -41.7, -31.8, -34.2, -38.6, -31.2, -35.3, -30.1, -28.3, -27.3, -35.0, -25.8, -19.6,
        -16.2, -6.0, -13.0, -3.1, 0.5, -7.0, 1.7]
errors = [5., 7., 5.0, 1.0, 5.0, 7.0, 1.0, 5.0, 1.0, 5.0, 6.0, 2.0,
        1.4, 4.0, 5.0, 6.0, 2.0, 5.0, 1.3,
        4.0, 4.6, 2.0, 4.0, 2.0, 1.3, 2.0, 2.0, 5.0, 5.0, 4.2, 4.2,
        5.0, 4.0, 3.0, 4.0, 4.0, 2.0, 4.0]

color1 = 'blue'
color2 = 'red'

for temperature1, temperature2, svc, error in zip(temperatures, temperatures[1:], svcs, errors):
    plt.scatter(temperature1, svc, marker = '*', color = 'k')

    if temperature1 == temperature2:
        plt.errorbar(temperature1, svc, error, color = color1)
    else:
        plt.errorbar(temperature1, svc, error, color = color2)
#plt.errorbar(213., -86.3, yerr = 5., color = 'blue')
#plt.errorbar(213., -94.0, yerr = 7., color = 'red')
plt.show()

#plt.errorbar(xrange(5), [2, 5, 3, 4, 7], yerr = [[1, 4, 2, 3, 6], [4, 10, 6, 8, 14]])
