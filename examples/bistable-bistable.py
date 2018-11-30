# Copyright (C) 2018 Urbain Vaes

import hermipy
from hermipy.equations import McKean_Vlasov as equation
import sympy as sym

# Configuration dicionaries
misc, eq, num = {}, {}, {}

# Variables and function
x, y, f = equation.x, equation.y, equation.f

# Short-hand notation
r = sym.Rational

# For bifurcation
# num['degree'] = 120  # degree of approximation
# num['n_points_num'] = 2*num['degree'] + 1  # (*2 for varf)
# num['μx'] = r(0, 4)
# num['μy'] = r(0, 4)
# num['σx'] = r(1, 50)
# num['σy'] = r(1, 10)
# num['λ'] = r(1, 2)
# num['index_set'] = 'rectangle'

# For convergence
num['degree'] = 250  # Or 250
num['n_points_num'] = 350
num['μx'] = r(0, 4)
num['μy'] = r(0, 4)
num['σx'] = r(1, 100)
num['σy'] = r(1, 10)
num['λ'] = r(1, 2)
num['index_set'] = 'rectangle'

# Scalar parameters of the equation
eq['ε'] = r(1, 4)
eq['γ'] = r(0)
eq['θ'] = r(1)

# For index set
# num['degree'] = 20  # degree of approximation
# num['n_points_num'] = 2*num['degree'] + 1  # (*2 for varf)
# num['μx'] = r(0, 1)
# num['μy'] = r(0, 1)
# num['σx'] = r(1, 15)
# num['σy'] = r(1, 15)
# num['λ'] = r(1, 2)
# num['index_set'] = 'triangle'
# eq['ε'] = r(1, 100)
# eq['γ'] = r(0)
# eq['θ'] = r(1)

# Functional parameters of the equation
eq['Vp'] = x**4/4 - x**2/2
eq['Vy'] = y**4/4 - y**2/2

# Miscellaneous parameters
misc['cache'] = True
misc['parallel'] = False
misc['tensorize'] = True
misc['sparse'] = True
misc['trails'] = False
misc['verbose'] = False
misc['symbolic'] = 0
misc['plots'] = False

# Asymptotic solution
# degree, Vy = 150, eq['Vy']
# hermipy.settings['cache'] = True
# fy = sym.Function('f')(y)
# qy = hermipy.Quad.gauss_hermite(2*degree+1, mean=[0], cov=[[.05]], dirs=[1])
# integral = qy.integrate(sym.exp(-Vy), flat=True)
# factor = sym.sqrt(qy.position.weight() / (sym.exp(-Vy) / integral))
# qy.factor = hermipy.Function(factor, dirs=[1])
# gen = -Vy.diff(y)*fy.diff(y) + fy.diff(y, y)
# L0 = qy.discretize_op(gen, degree=degree)
# ty = qy.transform('y', degree=degree)
# vy = qy.varf('y', degree=degree)
# A = ty*((-L0).solve(vy((-L0).solve(vy((-L0).solve(ty, remove0=True)), remove0=True)), remove0=True))
# B = ty*((-L0).solve((-L0).solve(ty, remove0=True), remove0=True))
# C = ty*(-L0).solve(ty)
# A, B, C = float(A), float(B), float(C)
# print(A/C**2, B/C)
