from sympy import symbols, sin, cos, diff

theta, varphi = symbols('theta varphi')

expr = sin(theta)**2 * cos(varphi)**2

z = sin(theta)
print diff(expr, z)