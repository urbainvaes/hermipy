import equation
import sympy as sym


def gaussian(mean, var):
    return r(.5)*(x - mean)*(x - mean)/(var)


# Configuration dicionaries
glob, params, functions, numerics = {}, {}, {}, {}
sym_params, sym_functions = {}, {}
equation.fill_params(sym_params, sym_functions)

# Short-hand notation
sp = sym_params

# Variables and function
x, y, f = equation.x, equation.y, equation.f

# Short-hand notation
r = sym.Rational

# Configuration of numerical method
numerics['degree'] = 50  # degree of approximation
numerics['n_points_num'] = 2*numerics['degree'] + 1  # (*2 for varf)

# Parameters of the equation
params['βx'] = r(1)
params['βy'] = r(1)
params['ε']  = r(.5)
params['γ']  = r(0)
params['θ']  = r(0)

# Parameters of the potential in the x equation
sym_params['mx'] = sym.symbols('mx', real=True)
sym_params['sx'] = sym.symbols('sx', real=True, positive=True)
params['mx'] = 0
params['sx'] = 1
functions['Vp'] = gaussian(sp['mx'], sp['sx'])/sp['βx']
functions['Vp'] = x**4/4 - x**2/2

# Parameters of approximating potential
sym_params['μx'] = sym.symbols('μx', real=True)
sym_params['σx'] = sym.symbols('σx', real=True, positive=True)
params['μx'] = r(1, 5)
params['σx'] = r(1, 10)
functions['Vq'] = gaussian(sp['μx'], sp['σx'])/sp['βx']

# Miscellaneous parameters
glob['cache'] = False
glob['symbolic'] = False
