#!/usr/bin/env python

# IMPORT MODULES {{{

import math
import sympy as sy
import sympy.printing as syp
import sympy.plotting as splot
import numpy as np
import numpy.linalg as la
import scipy.sparse
import scipy.sparse.linalg as las
import matplotlib.pyplot as plt

from libhermite import hermite as hm
sy.init_printing()

# }}}
# SYMBOLIC PARAMETERS OF STOCHASTIC SYSTEM {{{

# Variables
t, x, y = sy.symbols('t x y', real=True)

# Real parameters
beta_x, beta_y = sy.symbols('βx βy', real=True, positive=True)
epsilon = sy.symbols('ε', real=True, positive=True)
gamma = sy.symbols('γ', real=True, positive=True)
theta = sy.symbols('θ', real=True, positive=True)
m = sy.symbols('m', real=True, positive=True)

# Functional parameters
potential_p = sy.Function('V')(x)
potential_q = sy.Function('Vq')(x)
potential_y = sy.Function('Vy')(y)

# Dictionaries to store parameters
params, functions = {}, {}

# }}}
# DATA AND PARAMETERS FOR NUMERICAL SIMULATION {{{

# Degree of approximation
degree = 50

# Number of points in quadrature (*2 for varf)
n_points_num = 2*degree + 1

# Real parameters of the stochastic system
params[beta_x]  = sy.Rational(1)
params[theta]   = sy.Rational(0)
params[beta_y]  = sy.Rational(1)
params[epsilon] = sy.Rational(1, 2)
params[gamma]   = 0

# Potential
cov_p =  sy.symbols('sx')
# params[cov_p] = sy.Rational(1, 1)
# functions[potential_p] = sy.Rational(1, 2)*(x)*(x)/(beta_x*cov_p)
functions[potential_p] = x**4/4 - x**2/2
# potential_p = x**2/2 + 10*sy.cos(x)

# Parameters for approximating potential
mean_x, cov_x = sy.symbols('μx σx')
params[mean_x] = sy.Rational(1, 5)
params[cov_x] = sy.Rational(1, 10)
functions[potential_q] = sy.Rational(0.5)*(x - mean_x)*(x - mean_x)/(beta_x*cov_x)

# Potential of y process (Use scaling such that the coefficient of the asymptotic noise is 1)
alpha = sy.sqrt(2/beta_y)
cov_y = 1/(alpha*beta_y)
functions[potential_y] = y*y/(2*beta_y*cov_y)

# Mapping coefficient
map_quad = 2

# }}}
# ABSTRACT SYMBOLIC CALCULATIONS {{{

real_params = [(k, params[k]) for k in params]
functional_params = [(k, functions[k]) for k in functions]

# Solution of the Fokker-Planck equation
r = sy.Function('ρ')(x, y, t)

# Mapped solution: u = ρ * e^{V/2} * e^{Vq/2} to equation with BK operator
u = sy.Function('u')(x, y, t)
u_xy = sy.Function('u')(x, y)
u_x = sy.Function('u')(x)
u_y = sy.Function('u')(y)

# Fokker-Planck operator associated with potential
def forward(f):
    d = sy.diff
    drift_x = d(d(potential_p, x) * f + theta * (x - m) * f + (1 - gamma) * sy.sqrt(2 / beta_x) * y * f / epsilon, x)
    diff_x = gamma * (1/beta_x) * d(d(f, x), x)
    drift_y = (1/epsilon**2) * d(d(potential_y, y) * f, y)
    diff_y = (1/epsilon**2) * (1/beta_y) * d(d(f, y), y)
    return drift_x + diff_x + drift_y + diff_y

# Factors to map Fokker-Planck to backward Kolmogorov operator
factor_y = sy.exp(- beta_y * potential_y)
factor_q = sy.exp(- beta_x * potential_q)
factor_p = sy.exp(- beta_x * potential_p)

if map_quad == 0:
    factor_x = sy.exp(- beta_x * potential_p)
elif map_quad == 1:
    factor_x = sy.exp(- beta_x * potential_q)
elif map_quad == 2:
    factor_x = sy.exp(- beta_x/2 * (potential_q + potential_p))

factor = factor_x * factor_y

# Fokker-Planck and BK equations to solve (= 0)
fk = sy.diff(r, t) - forward(r)
bk = sy.simplify(fk.subs(r, u*factor).doit()/factor)
operator = - bk.subs(u, u_xy).doit()

# }}}
# {{{ CALCULATE THE EXACT SOLUTION

def solve():
    if functions[potential_p].diff(x, x, x) != 0:
        return False
    cx, cxy, cy = sy.symbols('cx cxy cy')
    def fk_expand(f):
        result = forward(f).subs(functional_params).doit()
        return result.subs(gamma, params[gamma])
    ansatz = sy.exp(-(cx*x*x + cxy*x*y + cy*y*y))
    res = (fk_expand(ansatz)/ansatz).doit().cancel()
    coeffs_dict = res.as_poly(x, y).as_dict()
    coeffs_list = list(coeffs_dict[k] for k in coeffs_dict)
    (sol_cx, sol_cxy, sol_cy) = sy.solve(coeffs_list, cx, cxy, cy)[0]
    solution = ansatz.subs([(cx, sol_cx), (cxy, sol_cxy), (cy, sol_cy)])
    assert fk_expand(solution).doit().expand() == 0
    return solution

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

cov_y = float(cov_y.subs(real_params))
factor_x = factor_x.subs(functional_params).subs(real_params)
factor_y = factor_y.subs(functional_params).subs(real_params)
factor   = factor.subs(functional_params).subs(real_params)
operator = operator.subs(functional_params).subs(real_params)
operator = operator.doit().expand()

for key in params:
    params[key] = float(params[key])

for key in functions:
    functions[key] = functions[key].subs(functional_params).subs(real_params)

solution = solve()
if solution:
    solution = solution.subs(real_params)
    splot.plot3d(solution.subs(real_params), (x, -1, 1), (y, -1, 1))

m_operator = operator.diff(m)
r_operator = (operator - m * m_operator).cancel()

# }}}
# DEFINE QUADRATURE FOR NUMERICS {{{

# For numerics
degrees = np.arange(degree + 1)
new_q = hm.Quad.gauss_hermite
mean_xy, cov_xy = [params[mean_x], 0], [[params[cov_x], 0], [0, cov_y]]
quad_num = new_q(n_points_num, dim=2, mean=mean_xy, cov=cov_xy)

# Calculate limits of resolution
band_width = np.sqrt(2) * np.sqrt(2*degree + 1)
x_min = params[mean_x] - band_width * np.sqrt(params[cov_x])
x_max = params[mean_x] + band_width * np.sqrt(params[cov_x])
y_min = 0 - band_width * np.sqrt(cov_y)
y_max = 0 + band_width * np.sqrt(cov_y)

# }}}
# DEFINE QUADRATURES FOR VISUALIZATION {{{

nv = 200
# bounds_x, bounds_y = band_width*np.sqrt(cov_x), band_width*np.sqrt(cov_y)
# bounds_x, bounds_y = 5*np.sqrt(cov_x), 5*np.sqrt(cov_y)
bounds_x, bounds_y = 3, 4

quad_num_x = new_q(n_points_num, dim=1, mean=[params[mean_x]], cov=[[params[cov_x]]])
quad_visu_xy = hm.Quad.newton_cotes([nv, nv], [bounds_x, bounds_y])
quad_visu_x = hm.Quad.newton_cotes([nv], [bounds_x])
quad_visu_y = hm.Quad.newton_cotes([nv], [bounds_y])

factor_visu_xy = quad_visu_xy.discretize(factor)
factor_visu_x = quad_visu_x.discretize(factor_x)
factor_visu_y = quad_visu_y.discretize(factor_y.subs(y, x))

x_visu = quad_visu_x.discretize('x')
y_visu = quad_visu_y.discretize('x')  # x is the default variable name

# }}}
# CALCULATE ASYMPTOTIC SOLUTION {{{
asymptotic_sol = sy.exp(- params[beta_x] * functions[potential_p])
series_asymptotic_sol = quad_num_x.transform(asymptotic_sol/factor_x, degree)
norm_asymptotic_sol = la.norm(series_asymptotic_sol.coeffs, 2)
series_asymptotic_sol.coeffs = series_asymptotic_sol.coeffs / norm_asymptotic_sol
asymptotic_solution_visu = quad_visu_x.eval(series_asymptotic_sol, degree)*factor_visu_x
disc_asymptotic_sol = quad_visu_x.discretize(asymptotic_sol/norm_asymptotic_sol)
# }}}
# SPECTRAL METHOD FOR STATIONARY EQUATION {{{

m_mat_operator =  quad_num.discretize_op(m_operator, u_xy, degree, 2)
r_mat_operator =  quad_num.discretize_op(r_operator, u_xy, degree, 2)

m_values = np.linspace(-1, 1, 11)

# m_num = 0
# for i in range(10):
#     total_mat = r_mat_operator + m_num * m_mat_operator


# Calculate eigenvector in kernel
print("Solving the eigenvalue problem...")
mat_operator = r_mat_operator
sparse_operator = scipy.sparse.csr_matrix(mat_operator)
asymptotic_sol_2d = sy.exp(- params[beta_x] * functions[potential_p]
                           - params[beta_y] * functions[potential_y])
v0 = quad_num.transform(asymptotic_sol_2d/factor, degree, norm=True).coeffs
eigen_values, eigen_vectors = las.eigs(mat_operator, v0=v0, k=3, which='LR', ncv=10)

print(eigen_values)
for eigen_value, eigen_vector in zip(eigen_values, eigen_vectors.T):
    solution = np.real(eigen_vector)

    solution_xy = solution * np.sign(solution[0])
    solution_x = hm.project(solution_xy, 2, 0)
    solution_y = hm.project(solution_xy, 2, 1)

    series_xy = hm.Series(solution_xy, mean=mean_xy, cov=cov_xy, norm=True)
    series_x = hm.Series(solution_x, mean=[params[mean_x]], cov=[[params[cov_x]]], norm=True)
    series_y = hm.Series(solution_y, mean=[0], cov=[[cov_y]], norm=True)

    solution_visu_xy = quad_visu_xy.eval(series_xy, degree)*factor_visu_xy
    solution_visu_xy = solution_visu_xy.reshape(nv, nv).T
    solution_visu_x = quad_visu_x.eval(series_x, degree) * factor_visu_x
    solution_visu_y = quad_visu_y.eval(series_y, degree) * factor_visu_y

    cont = plt.contourf(x_visu, y_visu, solution_visu_xy, 100)
    plt.title("Eigenvalue: " + str(eigen_value) + ", Minimum: " + str(np.min(solution_visu_xy)))
    plt.colorbar(cont)
    plt.show()

fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)

ax11.plot(x_visu, solution_visu_x, label="Coloured noise")
ax11.plot(x_visu, disc_asymptotic_sol, label="Asymptotic solution")
ax11.set_title("Comparison of the solutions")
ax11.legend()

# solution = solve()
# if solution:
#     ax12 = splot.plot3d(solution.subs(real_params), (x, -1, 1), (y, -1, 1))
# else:

ax12.bar(degrees, series_x.coeffs)
ax12.set_title("Hermite coefficients of numerical solution")

ax21.bar(degrees, series_x.coeffs - series_asymptotic_sol.coeffs)
ax21.set_title("Difference between Hermite coefficients")

cont = ax22.contourf(x_visu, y_visu, solution_visu_xy, 100)
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

# fig, ax = plt.subplots(1, 1)
# cont = ax.contourf(x_visu, y_visu, solutions_visu[0] - solutions_visu[1], 20)
# plt.colorbar(cont)
# plt.show()

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
