# Copyright (C) 2018 Urbain Vaes

from hermipy.equations import McKean_Vlasov as equation
import sympy as sym

# Configuration dicionaries
misc, eq, num = {}, {}, {}

# Variables and function
x, y, f = equation.x, equation.y, equation.f

# Short-hand notation
r = sym.Rational

# Configuration of numerical method
num['degree'] = 100  # degree of approximation
num['n_points_num'] = 2*num['degree'] + 1  # (*2 for varf)
num['μx'] = r(0, 4)
num['μy'] = r(0, 4)
num['σx'] = r(1, 20)
num['σy'] = r(1, 1)
num['λ'] = r(1, 2)
num['index_set'] = 'cube'

# Scalar parameters of the equation
# eq['β'] = r(1)
eq['θ'] = r(1)
eq['γ'] = r(0)

# Functional parameters of the equation
eq['Vp'] = x**4/4 - x**2/2
eq['Vy'] = y**2/2

# Miscellaneous parameters
misc['cache'] = True
misc['tensorize'] = True
misc['sparse'] = True
misc['trails'] = False
misc['verbose'] = False
misc['symbolic'] = 0  # Values 0, 1, 2
