from jitcode import jitcode, provide_basic_symbols
import numpy as np

a = -0.025794
b1 = 0.0065
b2 = 0.0135
c = 0.02
k = 0.128

t, y = provide_basic_symbols()
f = [
	y(0) * (a - y(0)) * (y(0) - 1.0) - y(1) + k * (y(2) - y(0)),
	b1 * y(0) - c * y(1),
	y(2) * (a - y(2)) * (y(2) - 1.0) - y(3) + k * (y(0) - y(2)),
	b2 * y(2) - c * y(3)
]

initial_state = np.array([1., 2., 3., 4.])

ODE = jitcode(f)
ODE.set_integrator("dopri5")
ODE.set_initial_value(initial_state, 0.0)

times = range(10, 100000, 10)
data = []
for time in times:
	data.append(ODE.integrate(time))

np.savetxt("timeseries.dat", data)



