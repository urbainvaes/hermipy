# IMPORT MODULES {{{
import sympy as sy
import numpy as np
import spectral as sp

# import os
# import sys
# import re
import time
# import sympy.printing
# import sympy.parsing.sympy_parser as symparser

import matplotlib.pyplot as plt

import matplotlib
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

# potential = x**4/4 - x**2/2
# potential = x**4
potential = x**2/2 + 10*sy.cos(x)

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

# }}}
# QUADRATIC POTENTIAL FOR APPROXIMATION {{{

mean = 0
cov = .1

potential_quad = 0.5*np.log(2*np.pi*cov) + (x - mean)*(x - mean)/(2 * cov)
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

# Factor between backward Kolmogorov and Schrődinger equation
sqrt_gaussian = sy.sqrt(rho_gaussian)
print("Linear term corresponding to quadratic potential:")
print(linear_gaussian)

# }}}
# ---- NUMERICAL METHOD ----
sy.init_printing()

# Number of degrees in Hermite expansion
degree = 30
degrees = np.arange(degree + 1)

# Fine quadrature for function evaluation
n_points_fine = 200
quad_fine = sp.Quad(n_points_fine, dim=1, mean=[mean], cov=[[cov]])
x_fine = quad_fine.discretize('x')
sqrt_gaussian_fine = quad_fine.discretize(sqrt_gaussian)

# Coarse quadrature for Hermite transform
n_points_coarse = degree + 1
quad_coarse = sp.Quad(n_points_coarse, dim=1, mean=[mean], cov=[[cov]])
sqrt_gaussian_coarse = quad_coarse.discretize(sqrt_gaussian)

# Difference between linear terms
diff_linear = linear - linear_gaussian

# Print approximating functions to Schrődinger equation
fig, ax = plt.subplots(1, 1)
for i in degrees:
    h = np.zeros(degree + 1)
    h[i] = 1
    Eh = quad_fine.eval(h, degree)[0] * sqrt_gaussian_fine
    ax.plot(x_fine, Eh)
plt.show()

# Discretize functions


# Time step and number of iterations
dt = 2e-2
n_iter = 100000

# Eigenvalues of the operator
# eigenvalues_gaussian = - np.arange(degree + 1)
eigenvalues_gaussian = np.zeros(degree + 1) - 1.

# Initial condition
u_init = 1
u_coarse = quad_coarse.discretize(u_init)
diff_linear_coarse = quad_coarse.discretize(diff_linear)
h_n = quad_coarse.transform(u_coarse, degree)

# Exact solution on fine grid
u_exact_fine = quad_fine.discretize(sy.sqrt(rho))
u_exact_coarse = quad_coarse.discretize(sy.sqrt(rho))
u_approx_fine = quad_fine.eval(quad_coarse.transform(u_exact_coarse/sqrt_gaussian_coarse, degree), degree)[0]*sqrt_gaussian_fine
norm_aux = quad_coarse.transform(u_exact_coarse/sqrt_gaussian_coarse, degree)
norm = np.sqrt(np.sum(np.square(norm_aux)))
u_exact_fine = u_exact_fine / norm
u_approx_fine = u_approx_fine / norm

plt.plot(x_fine, u_exact_fine)
plt.show()

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1)

# Activate interactive plotting
plt.ion()

for i in range(n_iter):

    # Plotting {{{
    if i % 1 == 0:
        plt.pause(.01)

        # Representation of u on fine grid
        u_fine = quad_fine.eval(h_n, degree)[0]

        # Plot solution in real space
        ax1.clear()
        ax1.set_title("Solution to Schrődinger equation")
        ax1.plot(x_fine, u_fine * sqrt_gaussian_fine)
        ax1.plot(x_fine, u_exact_fine)
        ax1.plot(x_fine, u_approx_fine)

        # Plot Hermite transform
        ax2.clear()
        ax2.set_title("Hermite coefficients of the solution")
        ax2.bar(degrees, h_n)

        plt.draw()
    # }}}

    # h_n = quad.transform(u_n, degree)
    h_n = h_n + dt/2 * eigenvalues_gaussian * h_n
    u_coarse = quad_coarse.eval(h_n, degree)[0]
    h_n = h_n + dt * quad_coarse.transform(u_coarse*diff_linear_coarse, degree)
    h_n = h_n + dt/2 * eigenvalues_gaussian * h_n

    # Normalization
    h_n = h_n/np.sqrt(np.sum(np.square(h_n)))

    # u_n = u_n + dt * diff_linear_n * u_n

# # plt.plot(points[0], potential_n)
