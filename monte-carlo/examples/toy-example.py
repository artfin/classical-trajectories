import math
import random

# out function to integrate
def f(x):
	return math.sin(x)

# define any xmin-xmax interval
xmin = 0.0
xmax = 2.0 * math.pi

# fin ymin-ymax
numSteps = 1000
ymin = f(xmin)
ymax = ymin
for i in xrange(numSteps):
	x = xmin + (xmax - xmin) * float(i) / numSteps
	y = f(x)
	if y < ymin: ymin = y
	if y > ymax: ymax = y

print 'ymin: {0}; ymax: {1}'.format(ymin, ymax)

# Monte Carlo
rectArea = (xmax - xmin) * (ymax - ymin)
numPoints = 100000000
ctr = 0

for j in xrange(numPoints):
	x = xmin + (xmax - xmin) * random.random()
	y = ymin + (ymax - ymin) * random.random()
	if math.fabs(y) <= math.fabs(f(x)):
		if f(x) > 0 and y > 0 and y <= f(x):
			ctr += 1 # area over x-axis is positive
		if f(x) < 0 and y <0 and y >= f(x):
			ctr -= 1 # area under x-axis is negative

fnArea = rectArea * float(ctr) / numPoints
print 'Numerical integration = ' + str(fnArea)



