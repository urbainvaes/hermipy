#!/usr/bin/env python

# IMPORT MODULES {{{

import itertools
import numpy.linalg as la
import scipy.sparse.linalg as las
import scipy.special
import scipy.sparse
import sympy as sy
import sympy.printing as syp
import numpy as np
import matplotlib.pyplot as plt
import math

from libhermite import hermite_python as hm
sy.init_printing()

# }}}
# PARAMETERS FOR NUMERICAL SIMULATION {{{

# Number of space dimensions
dim = 1

# Degree of approximation
degree = 20

# Number of points in quadrature (*2 for varf)
n_points_num = 2*degree + 1

# Parameters of the stochastic system
beta_x = 10
epsilon = 2**-1

# Parameters of the approximating Gaussian
mean_x = .2
cov_x = .1

# Potential
x = sy.symbols('x')
potential_p = x**4/4 - x**2/2
# potential_p = x*x/(2*beta_x*cov_x)
# potential_p = x**2/2 + 10*sy.cos(x)

# Mapping coefficient
map_quad = 2

# }}}
# ABSTRACT SYMBOLIC CALCULATIONS {{{

# Time, inverse temperatures, fraction of noise that is white, and parameter
# determining how close to white noise we are
t, beta_xa, epsilon_a = sy.symbols('t βx ε')

# Potential of interest
potential_pa = sy.Function('V')(x)

# Quadratic potential used for the approximation
potential_qa = sy.Function('Vq')(x)

# Solution of the Fokker-Planck equation
r = sy.Function('ρ')(x, t)

# Mapped solution: u = ρ * e^{V/2} * e^{Vq/2} to equation with BK operator
u = sy.Function('u')(x, t)
u_x = sy.Function('u')(x)


# Fokker-Planck operator associated with potential
def forward(potential, f):
    d = sy.diff
    drift_x = d(d(potential, x) * f, x)
    diff_x = (1/beta_xa) * d(d(f, x), x)
    return drift_x + diff_x


# Factors to map Fokker-Planck to backward Kolmogorov operator
factor_qa = sy.exp(- beta_xa * potential_qa)
factor_pa = sy.exp(- beta_xa * potential_pa)

if map_quad == 0:
    factor_a = sy.exp(- beta_xa * potential_pa)
elif map_quad == 1:
    factor_a = sy.exp(- beta_xa * potential_qa)
elif map_quad == 2:
    factor_a = sy.exp(- beta_xa/2 * (potential_qa + potential_pa))

# Fokker-Planck and BK equations to solve (= 0)
fk = sy.diff(r, t) - forward(potential_pa, r)
bk = sy.simplify(fk.subs(r, u*factor_a).doit()/factor_a)
operator_rhs_a = - bk.subs(u, u_x).doit()

# }}}
# PRINT TO STDOUT {{{

print_out = True

if print_out:

    print("Fokker-Planck equation to solve: ")
    syp.pprint(fk)

    print("Mapping to an equation with BK operator")
    syp.pprint(bk)

# }}}
# EVALUATE ABSTRACT EXPRESSIONS FOR PROBLEM AT HAND {{{

# Quadratic potential for approximation
potential_q = (x - mean_x)*(x - mean_x)/(2*beta_x*cov_x)


def evaluate(sym):
    sym = sym.subs(beta_xa, beta_x)
    sym = sym.subs(epsilon_a, epsilon)
    sym = sym.subs(potential_pa, potential_p)
    sym = sym.subs(potential_qa, potential_q)
    return sym


factor = evaluate(factor_a)
operator_rhs = evaluate(operator_rhs_a).doit()

# }}}
# SPLITTING OF THE OPERATOR {{{


def multi_indices(dim, deg_max, deg_min=0):
    return [m for m in itertools.product(range(deg_max+1), repeat=dim) if
            sum(m) <= deg_max and sum(m) >= deg_min]


def split_operator(op, func, variables):
    result, rem, order = [], op, 2
    for m in multi_indices(len(variables), order):
        test, der = 1, func
        for i, v in zip(m, variables):
            test *= v**i/math.factorial(i)
            der = sy.diff(der, v, i)
        term = rem.subs(func, test).doit()
        rem = (rem - term*der).simplify()
        result.append(term.simplify())
    # assert rem == 0
    return result


split_op = split_operator(operator_rhs, u_x, [x])

# }}}
# DEFINE QUADRATURE FOR NUMERICS {{{

# For numerics
degrees = np.arange(degree + 1)
new_q = hm.Quad.gauss_hermite
quad_num = new_q(n_points_num, dim=1, mean=[mean_x], cov=[[cov_x]])

# Calculate limits of resolution
band_width = np.sqrt(2) * np.sqrt(2*degree + 1)
x_min = mean_x - band_width * np.sqrt(cov_x)
x_max = mean_x + band_width * np.sqrt(cov_x)

# }}}
# SPECTRAL METHOD FOR STATIONARY EQUATION {{{

n_polys = int(scipy.special.binom(degree + dim, degree))
mat_operator = np.zeros((n_polys, n_polys))
mult = list(multi_indices(dim, 2))
for m, coeff in zip(mult, split_op):
    mat_operator += quad_num.varfd(coeff, degree, [x]*m[0])

# Calculate eigenvector in kernel
print("Solving the eigenvalue problem...")
sparse_operator = scipy.sparse.csr_matrix(mat_operator)
eigen_values, eigen_vectors = las.eigs(mat_operator, k=1, which='SM')
solution = np.real(eigen_vectors.T[0])
solution = solution / la.norm(solution, 2)
solution = solution * np.sign(solution[0])
series = hm.Series(solution, mean=[mean_x], cov=[[cov_x]])

# }}}
# QUADRATURES FOR VISUALIZATION {{{

nv = 200
bounds_x = 3
quad_visu = hm.Quad.newton_cotes([nv], [bounds_x])
factor_visu = quad_visu.discretize(factor)
x_visu = quad_visu.discretize('x')

# }}}
# PLOT OF THE SOLUTION OF STATIONARY FP {{{

exact_sol = sy.exp(- beta_x * (x**4/4 - x**2/2))
H_exact_sol = quad_num.transform(exact_sol/factor, degree)
norm_exact_sol = la.norm(H_exact_sol.coeffs, 2)
H_exact_sol.coeffs = H_exact_sol.coeffs / norm_exact_sol

approx_solution_visu = quad_visu.eval(series, degree)*factor_visu
exact_solution_visu = quad_visu.eval(H_exact_sol, degree)*factor_visu
actual_solution_visu = quad_visu.discretize(exact_sol/norm_exact_sol)

fig, ax11 = plt.subplots(1, 1)
ax11.plot(x_visu, approx_solution_visu)
ax11.plot(x_visu, exact_solution_visu)
ax11.plot(x_visu, actual_solution_visu)

if x_visu[0] < x_min:
    ax11.axvline(x=x_min)

if x_visu[-1] > x_max:
    ax11.axvline(x=x_max)

plt.show()
