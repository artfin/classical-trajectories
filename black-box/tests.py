def fcdiff(func, vals, argnum, step = 0.00001):
	_vals_p, _vals_m = vals[:], vals[:]
	_vals_p[argnum] = _vals_p[argnum] + step
	_vals_m[argnum] = _vals_m[argnum] - step
	return (func(*_vals_p) - func(*_vals_m)) / (2 * step)

print 'automatic dH/dq: {0}'.format(dham_dq(q, p, Jx0, Jy0, Jz0))
print 'numeric dH/dq: {0}'.format(fcdiff(func = hamiltonian, vals = [q, p, Jx0, Jy0, Jz0], argnum = 0))

print 'automatic dH/dp: {0}'.format(dham_dp(q, p, Jx0, Jy0, Jz0))
print 'numeric dH/dp: {0}'.format(fcdiff(func = hamiltonian, vals = [q, p, Jx0, Jy0, Jz0], argnum = 1))

print 'automatic dH/djx: {0}'.format(dham_djx(q, p, Jx0, Jy0, Jz0))
print 'numeric dH/djx: {0}'.format(fcdiff(func = hamiltonian, vals = [q, p, Jx0, Jy0, Jz0], argnum = 2))

print 'automatic dH/djy: {0}'.format(dham_djy(q, p, Jx0, Jy0, Jz0))
print 'numeric dH/djy: {0}'.format(fcdiff(func = hamiltonian, vals = [q, p, Jx0, Jy0, Jz0], argnum = 3))

print 'automatic dH/djz: {0}'.format(dham_djz(q, p, Jx0, Jy0, Jz0))
print 'numeric dH/djz: {0}'.format(fcdiff(func = hamiltonian, vals = [q, p, Jx0, Jy0, Jz0], argnum = 4))
