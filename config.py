import equation
import sympy as sym


def gaussian(mean, var):
    return r(.5)*(x - mean)*(x - mean)/(var)


# Configuration dicionaries
glob, eq, functions, num = {}, {}, {}, {}
sym_params = equation.forward_params()

# Short-hand notation
sp = sym_params

# Variables and function
x, y, f = equation.x, equation.y, equation.f

# Short-hand notation
r = sym.Rational

# Configuration of numerical method
num['degree'] = 40  # degree of approximation
num['n_points_num'] = 2*num['degree'] + 1  # (*2 for varf)
num['μx'] = r(1, 5)
num['σx'] = r(1, 10)
num['λ'] = r(1, 2)

# Scalar parameters of the equation
eq['βx'] = r(1)
eq['βy'] = r(1)
eq['ε'] = r(1)
eq['γ'] = r(0)
eq['θ'] = r(0)

# Functional parameters of the equation
eq['Vp'] = x**4/4 - x**2/2

# Parameters of the potential in the x equation
# sym_params['mx'] = sym.symbols('mx', real=True)
# sym_params['sx'] = sym.symbols('sx', real=True, positive=True)
# eq['mx'] = 0
# eq['sx'] = 1
# functions['Vp'] = gaussian(sp['mx'], sp['sx'])/sp['βx']

# Miscellaneous parameters
glob['cache'] = True
glob['symbolic'] = 2
