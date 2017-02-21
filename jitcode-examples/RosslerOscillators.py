from jitcode import jitcode, provide_basic_symbols
import numpy as np
import sympy

# parameters
N = 500
omega = np.random.uniform(0.9, 1.0, N)
a = 0.165
b = 0.2
c = 10.0
k = 0.01

# get adjacency matrix of a small-world network
A = small_world_network(
	number_of_nodes = N,
	nearest_neighbours = 20,
	rewiring_probability = 0.1
	)

# generate differential equations
t, y = provide_basic_symbols()

sum_z = sympy.Symbol('sum_z')
j = sympy.Symbol('j')
helpers = [( sum_z, sympy.Sum( y(3 * j + 2), (j, 0, N-1)) )]

def small_world_network(number_of_nodes, nearest_neighbours, rewiring_probability):
		n = number_of_nodes
		m = nearest_neighbours//2
		
		A = np.zeros( (n,n), dtype=bool )
		for i in range(n):
			for j in range(-m,m+1):
				A[i,(i+j)%n] = True
		
		# rewiring
		for i in range(n):
			for j in range(j):
				if A[i,j] and (np.random.random() < rewiring_probability):
					A[j,i] = A[i,j] = False
					while True:
						i_new,j_new = np.random.randint(0,n,2)
						if A[i_new,j_new] or i_new==j_new:
							continue
						else:
							A[j_new,i_new] = A[i_new,j_new] = True
							break
		
		return A

def f():
	for i in range(N):
		coupling_term = sympy.Mul(
			k,
			sum( (y(3*j)-y(3*i)) for j in range(N) if A[i, j] ),
			evaluate = False
		)
		yield -omega[i] * y(3*i+1) - y(3*i+2) + coupling_term
		yield omega[i] * y(3*i) + a*y(3*i+1)
		coupling_term_2 = k * (sum_z-N*y(3*i+2))
		yield b + y(3*i+2) * (y(3*i) - c) + coupling_term_2

# integrate

initial_state = np.random.random(3*N)

ODE = jitcode(f, helpers = helpers, n = 3*N)
ODE.generate_f_C(simplify = False, do_cse = False, chunk_size = 150)
ODE.set_integrator('dopri5')
ODE.set_initial_value(initial_state, 0.0)

# data structure: x[0], v[0], z[0], x[1], ..., x[N], v[N], z[N]
data = np.vstack(ODE.integrate(T) for T in range(10, 100000, 10))
















