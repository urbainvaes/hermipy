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

from libhermite import hermite as hm
sy.init_printing()

# }}}
# PARAMETERS FOR NUMERICAL SIMULATION {{{

# Degree of approximation
degree = 15

# Number of points in quadrature (*2 for varf)
n_points_num = 2*degree + 1

# Parameters of the stochastic system
beta_x = 3.
beta_y = 2.
epsilon = 2**0

# noise = (1-gamma) * coloured + gamma * white
gamma = 0

# Parameters of the approximating Gaussian
mean_x = .2
cov_x = .3

# Potential
x, y = sy.symbols('x y')
potential_p = x**4/4 - x**2/2
# potential_p = x*x/(2*beta_x*cov_x)
# potential_p = x**2/2 + 10*sy.cos(x)

# Mapping coefficient
map_quad = 2

# }}}
# ABSTRACT SYMBOLIC CALCULATIONS {{{

# Number of space dimensions
dim = 2

# Time, inverse temperatures, fraction of noise that is white, and parameter
# determining how close to white noise we are
t, beta_xa, beta_ya, gamma_a, epsilon_a = sy.symbols('t βx βy γ ε', real=True, positive=True)

# Use correct scaling
alpha_a = sy.sqrt(2/beta_ya)

# Potential of interest
potential_pa = sy.Function('V')(x)

# Quadratic potential used for the approximation
potential_qa = sy.Function('Vq')(x)

# Potential for coloured noise
cov_ya = 1/(alpha_a*beta_ya)
potential_ya = sy.Function('Vy')(y)

# Solution of the Fokker-Planck equation
r = sy.Function('ρ')(x, y, t)

# Mapped solution: u = ρ * e^{V/2} * e^{Vq/2} to equation with BK operator
u = sy.Function('u')(x, y, t)
u_xy = sy.Function('u')(x, y)
u_x = sy.Function('u')(x)
u_y = sy.Function('u')(y)


# Fokker-Planck operator associated with potential
def forward(potential, f):
    d = sy.diff
    drift_x = d(d(potential, x) * f + (1 - gamma_a) * sy.sqrt(2 / beta_xa) * y * f / epsilon_a, x)
    diff_x = gamma_a * (1/beta_xa) * d(d(f, x), x)
    drift_y = (1/epsilon_a**2) * d(d(potential_ya, y) * f, y)
    diff_y = (1/epsilon_a**2) * (1/beta_ya) * d(d(f, y), y)
    return drift_x + diff_x + drift_y + diff_y


# Factors to map Fokker-Planck to backward Kolmogorov operator
factor_ya = sy.exp(- beta_ya * potential_ya)
factor_qa = sy.exp(- beta_xa * potential_qa)
factor_pa = sy.exp(- beta_xa * potential_pa)

if map_quad == 0:
    factor_xa = sy.exp(- beta_xa * potential_pa)
elif map_quad == 1:
    factor_xa = sy.exp(- beta_xa * potential_qa)
elif map_quad == 2:
    factor_xa = sy.exp(- beta_xa/2 * (potential_qa + potential_pa))

factor_a = factor_xa * factor_ya

# Fokker-Planck and BK equations to solve (= 0)
fk = sy.diff(r, t) - forward(potential_pa, r)
bk = sy.simplify(fk.subs(r, u*factor_a).doit()/factor_a)
operator_rhs_a = - bk.subs(u, u_xy).doit()

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

# Coefficient of the y drift
alpha = alpha_a.subs(beta_xa, beta_x).subs(beta_ya, beta_y)

# Quadratic potential for approximation
potential_q = (x - mean_x)*(x - mean_x)/(2*beta_x*cov_x)

# Potential of y process
potential_y = y*y/(2*beta_y*cov_ya)

def evaluate(sym):
    sym = sym.subs(potential_pa, potential_p)
    sym = sym.subs(potential_qa, potential_q)
    sym = sym.subs(potential_ya, potential_y)
    sym = sym.subs(beta_xa, beta_x)
    sym = sym.subs(beta_ya, beta_y)
    sym = sym.subs(gamma_a, gamma)
    sym = sym.subs(epsilon_a, epsilon)
    return sym


cov_y = float(evaluate(cov_ya))

factor_x = evaluate(factor_xa)
factor_y = evaluate(factor_ya)
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


split_op = split_operator(operator_rhs, u_xy, (x, y))

# }}}
# DEFINE QUADRATURE FOR NUMERICS {{{

# For numerics
degrees = np.arange(degree + 1)
new_q = hm.Quad.gauss_hermite
mean_xy, cov_xy = [mean_x, 0], [[cov_x, 0], [0, cov_y]]
quad_num = new_q(n_points_num, dim=2, mean=mean_xy, cov=cov_xy)

# Calculate limits of resolution
band_width = np.sqrt(2) * np.sqrt(2*degree + 1)
x_min = mean_x - band_width * np.sqrt(cov_x)
x_max = mean_x + band_width * np.sqrt(cov_x)
y_min = 0 - band_width * np.sqrt(cov_y)
y_max = 0 + band_width * np.sqrt(cov_y)

# }}}
# SPECTRAL METHOD FOR STATIONARY EQUATION {{{

n_polys = int(scipy.special.binom(degree + dim, degree))
mat_operator = np.zeros((n_polys, n_polys))
mult = list(multi_indices(dim, 2))
for m, coeff in zip(mult, split_op):
    mat_operator += quad_num.dvarf(coeff, degree, [x]*m[0] + ['y']*m[1])

# Calculate eigenvector in kernel
print("Solving the eigenvalue problem...")
sparse_operator = scipy.sparse.csr_matrix(mat_operator)
eigen_values, eigen_vectors = las.eigs(mat_operator, k=1, which='SM')
solution = np.real(eigen_vectors.T[-1])

solution_xy = solution * np.sign(solution[0])
solution_x = hm.project(solution_xy, 2, 0)
solution_y = hm.project(solution_xy, 2, 1)

solution_xy = solution_xy / la.norm(solution_xy, 2)
solution_x = solution_x / la.norm(solution_x, 2)
solution_y = solution_y / la.norm(solution_y, 2)

series_xy = hm.Series(solution_xy, mean=mean_xy, cov=cov_xy)
series_x = hm.Series(solution_x, mean=[mean_x], cov=[[cov_x]])
series_y = hm.Series(solution_y, mean=[0], cov=[[cov_y]])

# assert(la.norm(series_y.coeffs[1:], 2) < 1e-2)

# }}}
# QUADRATURES FOR VISUALIZATION {{{

nv = 200
# bounds_x, bounds_y = band_width*np.sqrt(cov_x), band_width*np.sqrt(cov_y)
# bounds_x, bounds_y = 5*np.sqrt(cov_x), 5*np.sqrt(cov_y)
bounds_x, bounds_y = 3, 4

quad_num_x = new_q(n_points_num, dim=1, mean=[mean_x], cov=[[cov_x]])
quad_visu_xy = hm.Quad.newton_cotes([nv, nv], [bounds_x, bounds_y])
quad_visu_x = hm.Quad.newton_cotes([nv], [bounds_x])
quad_visu_y = hm.Quad.newton_cotes([nv], [bounds_y])

factor_visu_xy = quad_visu_xy.discretize(factor)
factor_visu_x = quad_visu_x.discretize(factor_x)
factor_visu_y = quad_visu_y.discretize(factor_y.subs(y, x))

x_visu = quad_visu_x.discretize('x')
y_visu = quad_visu_y.discretize('x')  # x is the default variable name

# }}}
# PLOT OF THE SOLUTION OF STATIONARY FP {{{

white_sol = sy.exp(- beta_x * potential_p)
series_white_sol = quad_num_x.transform(white_sol/factor_x, degree)
norm_white_sol = la.norm(series_white_sol.coeffs, 2)
series_white_sol.coeffs = series_white_sol.coeffs / norm_white_sol
white_solution_visu = quad_visu_x.eval(series_white_sol, degree)*factor_visu_x
disc_white_sol = quad_visu_x.discretize(white_sol/norm_white_sol)

solution_visu_xy = quad_visu_xy.eval(series_xy, degree)*factor_visu_xy
solution_visu_xy = solution_visu_xy.reshape(nv, nv).T
solution_visu_x = quad_visu_x.eval(series_x, degree) * factor_visu_x
solution_visu_y = quad_visu_y.eval(series_y, degree) * factor_visu_y

fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)

ax11.plot(x_visu, solution_visu_x, label="Coloured noise")
ax11.plot(x_visu, disc_white_sol, label="White noise")
ax11.set_title("Comparison of the solutions")
ax11.legend()

ax12.bar(degrees, series_x.coeffs)
ax12.set_title("Hermite coefficients of numerical solution")

ax21.bar(degrees, series_x.coeffs - series_white_sol.coeffs)
ax21.set_title("Difference between Hermite coefficients")

# ax21.plot(y_visu, solution_visu_y)
# ax21.set_title("Probability density of noise")

cont = ax22.contourf(x_visu, y_visu, solution_visu_xy, 20)
ax22.set_title("SM eigenvalues: " + str(eigen_values))

plt.colorbar(cont)

if x_visu[0] < x_min:
    ax11.axvline(x=x_min)
    ax22.axvline(x=x_min)

if x_visu[-1] > x_max:
    ax11.axvline(x=x_max)
    ax22.axvline(x=x_max)

if y_visu[0] < y_min:
    ax22.axhline(y=y_min)
    ax21.axvline(x=y_min)

if y_visu[-1] > y_max:
    ax21.axvline(x=y_max)
    ax22.axhline(y=y_max)

plt.show()

# }}}
# PLOT HERMITE FUNCTION OF HIGHEST DEGREE {{{

# fig, ax = plt.subplots(1, 1)
# ax.set_title("Hermite function of degree " + str(degree))
# ax.axvline(x=x_min)
# ax.axvline(x=x_max)
# ax.set_ylim((-2, 2))
# h_i = np.zeros(degree + 1)
# h_i[degree] = 1
# h_i = hm.Series(h_i, mean=[mean_x], cov=[[cov_x]])
# # Eh = quad_visu.eval(h_i, degree) * np.sqrt(factor_visu)
# Eh = quad_visu_x.eval(h_i, degree) * factor_visu_x
# ax.plot(x_visu, Eh)
# plt.show()

# }}}
# SOLUTION OF TIME-DEPENDENT EQUATION {{{

# # Time step and number of iterations
# dt = 2e-4*cov
# n_iter = 10**6

# # Initial condition
# u_init = factor_q / factor_p
# u_num = quad_num.discretize(u_init)
# Hu = quad_num.transform(u_num, degree)

# # Create figure with two subplots
# fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)

# # Activate interactive plotting
# plt.ion()

# # Splitting method
# for i in range(n_iter):

#     # Plotting {{{
#     if i % 1000 == 0:
#         plt.pause(.01)

#         # Representation of u on fine grid
#         u_visu = quad_visu.eval(Hu, degree)

#         # Plot solution in flat space
#         ax11.clear()
#         ax11.set_title("Solution to Schrődinger equation")
#         ax11.plot(x_visu, u_visu * factor_q_visu)
#         ax11.plot(x_visu, u_sol_visu * factor_q_visu)

#         # Plot solution in real space
#         ax12.clear()
#         ax12.set_title("Solution to Fokker-Planck equation")
#         ax12.plot(x_visu, u_visu * factor_visu)

#         # Plot Hermite transform
#         ax21.clear()
#         ax21.set_title("Hermite coefficients of the numerical solution")
#         ax21.bar(degrees, Hu.coeffs)

#         # Plot error for Hermite coefficients
#         ax22.clear()
#         ax22.set_title("Hermite coefficients of the stationary solution")
#         ax22.bar(degrees, Hu_sol.coeffs[0:degree + 1])

#         plt.draw()
#     # }}}

#     Hu.coeffs = Hu.coeffs + dt/2 * eigenvalues_gaussian * Hu.coeffs
#     u_num = quad_num.eval(Hu, degree)
#     Hu.coeffs = Hu.coeffs + dt * quad_num.transform(u_num*multiplication_num, degree).coeffs
#     Hu.coeffs = Hu.coeffs + dt/2 * eigenvalues_gaussian * Hu.coeffs

#     # Normalization
#     Hu.coeffs = Hu.coeffs/la.norm(Hu.coeffs, 2)

# import importlib
# importlib.reload(hm)
# }}}
