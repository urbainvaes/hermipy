# Copyright (C) 2018 Urbain Vaes

from hermipy.equations import McKean_Vlasov as equation
import sympy as sym
import hermipy

# Configuration dicionaries
misc, eq, num = {}, {}, {}

# Variables and function
x, y, f = equation.x, equation.y, equation.f

# Short-hand notation
r = sym.Rational

# Configuration of numerical method
num['degree'] = 80  # degree of approximation
num['n_points_num'] = 2*num['degree'] + 1  # (*2 for varf)
num['μx'] = r(0, 4)
num['μy'] = r(0, 4)
num['σx'] = r(1, 40)
num['σy'] = r(1, 10)
num['λ'] = r(1, 2)
num['index_set'] = 'rectangle'

# Scalar parameters of the equation
eq['ε'] = r(1, 20)
eq['γ'] = r(0)
eq['θ'] = r(1)

# Functional parameters of the equation
eq['Vp'] = x**4/4 - x**2/2

# Mean-zero
Z, m = 6.301119049538182, 0.8852269357209047
m = m + 0.152235
m = 0
m = r(m).limit_denominator(1e16)
eq['Vy'] = (y-m)**4/4 - (y-m)**2/2 + (y-m)

eq['Vy'] = y**4/4 - y**2/2
# eq['Vy'] = y**2/2

Vy = eq['Vy']
ny, μy, σy = num['n_points_num'], [num['μy']], [[num['σy']]]
qy = hermipy.Quad.gauss_hermite(ny, mean=μy, cov=σy, dirs=[1])
factor = sym.sqrt(qy.position.weight() * sym.exp(-Vy))
qy.factor = hermipy.Function(factor, dirs=[1])

fy = sym.Function('f')(y)
index_set, degree = num['index_set'], num['degree']
gen = (Vy.diff(y)*fy).diff(y) + fy.diff(y, y)

qy.factor = hermipy.Function(factor, dirs=[1])
L0 = qy.discretize_op(gen, degree=degree, index_set=index_set)
l, [e] = L0.eigs(k=1, which='LR')
vy = qy.varf('y', degree=degree, index_set=index_set)
coeff_noise = 1/sym.sqrt(sym.Rational(37243868, 52597017))

drift = - r((vy(e)*e).coeffs[0]).limit_denominator(1e16) * coeff_noise
num['drift_correction'] = drift

# import matplotlib.pyplot as plt
# import matplotlib
# import numpy as np
# matplotlib.rc('font', size=14)
# matplotlib.rc('font', family='serif')
# matplotlib.rc('text', usetex=True)
# def plot_log(x, y, file, lin=True):
#     x, y = np.asarray(x), np.asarray(y)
#     fig, ax = plt.subplots(1, 1)
#     ax.semilogy(x, y, 'k.')
#     ax.set_yscale('log', basey=2)
#     cut_off = 30
#     x_poly = x[0:cut_off + 1] if len(x) > cut_off else x
#     y_poly = y[0:cut_off + 1] if len(y) > cut_off else y
#     coeffs = np.polyfit(x_poly, np.log10(y_poly), 1)
#     if lin:
#         ax.semilogy(x, 10**coeffs[1] * 10**(coeffs[0]*x), 'r-',
#                     label='$y = C \\times 10^{{ {:.2f} \\, x}}$'.
#                           format(coeffs[0]))
#     plt.legend(loc='upper right')
#     plt.savefig(file, bbox_inches='tight')
#     plt.show()
# plot_log(degrees, np.abs(bias), "bias.eps")

# Z = qy.integrate(sym.exp(Vy), flat=True)
# m1 = qy.integrate(y*sym.exp(Vy), flat=True)
# assert abs(qy.integrate(y*sym.exp(Vy), flat=True)) < 1e-10

# qx.norm(qx.transform(rx, degree=d).eval(), qx.discretize(yx))
# qy.plot(sym.exp(-eq['Vy']))

# Miscellaneous parameters
misc['cache'] = True
misc['parallel'] = False
misc['tensorize'] = True
misc['sparse'] = True
misc['trails'] = False
misc['verbose'] = False
misc['symbolic'] = 0  # Values 0, 1, 2
misc['plots'] = False
