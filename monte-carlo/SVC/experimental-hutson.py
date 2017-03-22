import matplotlib.pyplot as plt

temperatures = [213., 223., 242., 262., 276., 288.2, 296., 303.2, 313.2, 323.1, 333.2, 365.]

split_rep = [-88.3, -80.2, -67.3, -56.3, -49.7, -44.7, -41.7, -39.1, -35.8, -32.8, -29.9, -22.0]

experiment = [-93.9, -84.5, -67.2, -53.8, -48.7, -40.3, -37.05, -34.21, -31.20, -28.30, -25.8, -19.6]

errors = [10., 10., 10., 10., 7., 7., 7., 7., 7., 7., 7., 7.]

for temperature, svc_hutson, svc, error in zip(temperatures, split_rep, experiment, errors):
    plt.scatter(temperature, svc, marker = '*', color = 'k')
    plt.scatter(temperature, svc_hutson, marker = 'x', color = 'magenta')
    plt.errorbar(temperature, svc, error, color = 'red')

#plt.show()
plt.savefig('hutson-article.png')

