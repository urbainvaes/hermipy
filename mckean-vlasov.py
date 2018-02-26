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

# Number of space dimensions
dim = 1

# Space variable
x = sy.symbols('x')

# Potential
potential = x**4/4 - x**2/2

# Number of discretization points
n_points = 100

# Nodes and weights of the Gauss-Hermite quadrature
points, weights = sp.hermegauss_nd(n_points, dim=dim)

# Initial condition
rho_init = 1/(2*sy.pi) * sy.exp(-x**2/2)

# Discretized initial condition
rho_init_n = sp.discretize(rho_init, points)

plt.plot(points[0], rho_init_n)
plt.show()


# import importlib
# importlib.reload(sp)
