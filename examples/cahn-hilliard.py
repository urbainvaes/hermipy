import hermipy.quad as hm
import sympy as sym

# Variables
x, y = sym.symbols('x y', real=True)

# Unknown function
f = sym.Function('f')(x, y)

# Initial condition
f_init = 2*(x*x + 4*y*y < 1) - 1

# Degree and number of points
degree = 50

# (*2 for varf)
n_points = 2*degree + 1

# Quadrature
quad = hm.Quad.gauss_hermite(n_points, dim=2)

# Heat operator
# operator = 
