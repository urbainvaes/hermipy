# IMPORT MODULES {{{
import sympy as sy
import numpy as np
import spectral as sp

import os
import sys
import re
import sympy
import sympy.printing
import sympy.parsing.sympy_parser as symparser

import matplotlib
import matplotlib.pyplot as plt
# }}}

import importlib
importlib.reload(sp)

# PARAMETERS OF THE EQUATION {{{

# Number of space dimensions
dim = 1

# Space variable
x = sy.symbols('x')

# Inverse temperature
beta = 1

potential = x**4/4 - x**2/2
# potential = x**4
# potential = x**2/2 + sy.cos(x)

dx_potential = sy.diff(potential, x)
dxx_potential = sy.diff(dx_potential, x)
rho = sy.exp(-beta*potential)


# Backward Kolmogorov operator
def backward(f):
    dx_f = sy.diff(f, x)
    dxx_f = sy.diff(dx_f, x)
    return - dx_potential * dx_f + (1/beta) * dxx_f


# Fokker-Planck operator
def forward(f):
    dx_f = sy.diff(f, x)
    dxx_f = sy.diff(dx_f, x)
    return - (dx_potential * dx_f + dxx_potential * f) * f + beta * dxx_f


# Linear term in Schrodinger equation
linear = sy.sqrt(rho) * backward(1/sy.sqrt(rho))
linear = sy.expand(sy.simplify(linear))

# Factor to pass from between Fokker-Planck to Schrodinger
factor = sy.exp(potential/2)
inv_factor = sy.exp(-potential/2)

# }}}
# QUADRATIC POTENTIAL FOR APPROXIMATION {{{

mean = 0
cov = 1

potential_quad = 0.5*sy.log(2*sy.pi*cov) * (x - mean)*(x - mean)/(2 * cov)
dx_potential_quad = sy.diff(potential_quad, x)
dxx_potential_quad = sy.diff(dx_potential_quad, x)
rho_gaussian = sy.exp(-beta*potential_quad)


# Backward Kolmogorov operator
def backward_quad(f):
    dx_f = sy.diff(f, x)
    dxx_f = sy.diff(dx_f, x)
    return - dx_potential_quad * dx_f + (1/beta) * dxx_f

# Linear term in the case of the quadratic potential
linear_gaussian = sy.sqrt(rho_gaussian) * backward_quad(1/sy.sqrt(rho_gaussian))
linear_gaussian = sy.expand(sy.simplify(linear_gaussian))

# }}}
# ---- NUMERICAL METHOD ----

# Difference between linear terms
diff_linear = linear - linear_gaussian

# Number of discretization points
n_points = 100

# Nodes and weights of the Gauss-Hermite quadrature
# points, weights = sp.hermegauss_nd(n_points, dim=dim)

u_init = 1

# Discretize functions
quad = sp.Quad(n_points, dim=1, mean=[mean], cov=[[cov]])
factor_n = quad.discretize(factor)
inv_factor_n = quad.discretize(inv_factor)
u_init_n = quad.discretize(u_init)
diff_linear_n = quad.discretize(diff_linear)

n_iter = 10
degree = n_points - 1

# Real x for plots
x_points = quad.discretize('x')
degrees = np.arange(degree + 1)


# Plot difference between linear terms
# plot = plt.plot(x_points, diff_linear_n)
# plt.show()


degree = 7
dt = 0.0
n_iter = 10
# Number of iterations

eigenvalues_gaussian = np.arange(degree + 1)
u_n = u_init_n
for i in range(n_iter):
    plot = plt.plot(x_points, u_n);
    plt.show(block=True)
    h_n = quad.transform(u_n, degree)
    h_n = h_n + dt * eigenvalues_gaussian * h_n
    u_n = quad.eval(h_n, degree)[0]
    u_n = u_n + dt * diff_linear_n * u_n

plt.show()

# # plt.plot(points[0], potential_n)
