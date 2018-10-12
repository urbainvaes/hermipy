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
num['degree'] = 100  # degree of approximation
num['n_points_num'] = 2*num['degree'] + 1  # (*2 for varf)
num['μx'] = r(0, 4)
num['μy'] = r(0, 4)
num['σx'] = r(1, 50)
num['σy'] = r(1, 10)
num['λ'] = r(1, 2)
num['index_set'] = 'cube'

# Scalar parameters of the equation
# eq['β'] = r(2**5, 2**2)
# eq['ε'] = r(1, 2**4)
# eq['β'] = r(2**6, 2**2)
# eq['β'] = r(2)
eq['ε'] = r(1, 20)
eq['γ'] = r(0)
eq['θ'] = r(1)
# eq['m'] = r(0)

# Functional parameters of the equation
eq['Vp'] = x**4/4 - x**2/2
# eq['Vp'] = x**4/4
# + (x - 1)**2/2
# eq['Vp'] = x**2/2
# eq['Vp'] = x**4

eq['Vy'] = y**4/4 - y**2/2
# eq['Vy'] = y**2/2


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
