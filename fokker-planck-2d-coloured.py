# IMPORT MODULES {{{

import numpy.linalg as la
import scipy.sparse.linalg as las
import sympy as sy
import sympy.printing as syp
import numpy as np
import matplotlib.pyplot as plt
import math

from libhermite import hermite_python as hm
sy.init_printing()

# }}}
# ABSTRACT SYMBOLIC CALCULATIONS {{{

# Number of space dimensions
dim = 2

# Space variable, time variable and inverse temperature
x, y, t, beta_xa, beta_ya, gamma_a, epsilon_a = sy.symbols('x y t βx βy γ ε')

# Use correct scaling
alpha_a = sy.sqrt(beta_xa/beta_ya)

# Potential of interest
potential_pa = sy.Function('V')(x)

# Quadratic potential used for the approximation
potential_qa = sy.Function('Vq')(x)

# Potential for coloured noise
cov_ya = 1/(alpha_a*beta_ya)
potential_ya = y*y/(2*beta_ya*cov_ya)

# Solution of the Fokker-Planck equation
r = sy.Function('ρ')(x, y, t)

# Mapped solution: u = ρ * e^{V/2} * e^{Vq/2} to equation with BK operator
u = sy.Function('u')(x, y, t)
u_xy = sy.Function('u')(x, y)
u_x = sy.Function('u')(x)
u_y = sy.Function('u')(y)

# Mapped solution of Schrodinger equation
v = sy.Function('v')(x, y, t)


# Fokker-Planck operator associated with potential
def forward(potential, f):
    d = sy.diff
    drift_x = d(d(potential, x) * f, x) + (1 - gamma_a) * y * f / epsilon_a
    diff_x = gamma_a * (1/beta_xa) * d(d(f, x), x)
    drift_y = (1/epsilon_a**2) * d(d(potential_ya, y) * f, y)
    diff_y = (1/epsilon_a**2) * (1/beta_ya) * d(d(f, y), y)
    return drift_x + diff_x + drift_y + diff_y


# Factors to map Fokker-Planck to backward Kolmogorov operator
factor_pa = sy.exp(- beta_ya * potential_ya/2) * \
            sy.exp(- beta_xa * potential_pa / 2)
factor_qa = sy.exp(- beta_ya * potential_ya/2) * \
            sy.exp(- beta_xa * potential_qa / 2)

# Mapping to L² space weighted by Gaussian
factor_a = factor_pa * factor_qa

# Fokker-Planck equation to solve (= 0)
fk = sy.diff(r, t) - forward(potential_pa, r)

# Backward Kolmogorov-like equation
bk = sy.simplify(fk.subs(r, u*factor_a).doit()/factor_a)

# Obtain RHS (if the left-hand side is the time derivative)
operator_rhs = - bk.subs(u, u_xy).doit()

# Obtain multiplication operator
multiplication_a = operator_rhs.subs(u_xy, 1).doit()
multiplication_ax, multiplication_ay = 0, 0

for term in multiplication_a.args:
    if x in term.free_symbols:
        multiplication_ax += term
    elif y in term.free_symbols:
        multiplication_ay += term

# Remainder
operator_xa = operator_rhs.subs(u_xy, u_x).doit() - multiplication_ay * u_x
operator_ya = operator_rhs.subs(u_xy, u_y).doit() - multiplication_ax * u_y
operator_xa, operator_ya = operator_xa.simplify(), operator_ya.simplify()
operator_xa, operator_ya = operator_xa.simplify(), operator_ya.simplify()

remainder = sy.simplify(operator_rhs
                        - operator_xa.subs(u_x, u_xy)
                        - operator_ya.subs(u_y, u_xy))

assert remainder == 0

# }}}
# PRINT TO STDOUT {{{
print("Fokker-Planck equation to solve: ")
syp.pprint(fk)

print("Mapping to an equation with BK operator")
syp.pprint(bk)

print("Operator in the right-hand side")
syp.pprint(operator_rhs)

print("X part of the operator")
syp.pprint(operator_xa)

print("Y part of the operator")
syp.pprint(operator_ya)
# }}}
# EVALUATE ABSTRACT EXPRESSIONS FOR PROBLEM AT HAND {{{

beta_x = 2.
beta_y = 0.5
epsilon = 0.01
gamma = 1

# Coefficient of the y drift
alpha = alpha_a.subs(beta_xa, beta_x).subs(beta_ya, beta_y)

# For numerical approximation
mean_x = .2
cov_x = .5

# potential_p = x**2/2 + 10*sy.cos(x)
potential_p = x**4/4 - x**2/2
potential_q = (x - mean_x)*(x - mean_x)/(2*beta_x*cov_x)


def evaluate(sym):
    sym = sym.subs(beta_xa, beta_x)
    sym = sym.subs(beta_ya, beta_y)
    sym = sym.subs(gamma_a, gamma)
    sym = sym.subs(epsilon_a, epsilon)
    sym = sym.subs(potential_pa, potential_p)
    sym = sym.subs(potential_qa, potential_q)
    return sym


cov_y = evaluate(cov_ya)

factor_p = evaluate(factor_pa)
factor_q = evaluate(factor_qa)
factor = factor_q*factor_p

factor_x, factor_y = 1, 1
for term in factor.args:
    if x in term.free_symbols:
        factor_x *= term
    elif y in term.free_symbols:
        factor_y *= term

operator_x = evaluate(operator_xa).doit().expand()
operator_y = evaluate(operator_ya).doit().expand()


def split_operator(op, func, var):
    result, rem, order = [], op, 2
    for i in range(order + 1):
        term = rem.subs(func, var**i/math.factorial(i)).doit()
        rem = (rem - term*sy.diff(func, var, i)).simplify()
        result.append(term.simplify())
    return result


split_operator_x = split_operator(operator_x, u_x, x)
split_operator_y = split_operator(operator_y, u_y, y)

split_operator_y = [term.subs(y, x) for term in split_operator_y]
factor_y = factor_y.subs(y, x)

# }}}
# DISCRETIZE VARIOUS FUNCTIONS ON GRID {{{

# For numerics
degree = 10
degrees = np.arange(degree + 1)
n_points_num = degree + 1
new_q = hm.Quad.gauss_hermite
quad_num_x = new_q(n_points_num, dim=1, mean=[mean_x], cov=[[cov_x]])
quad_num_y = new_q(n_points_num, dim=1, mean=[0], cov=[[cov_y]])

# Calculate limits of resolution
band_width = np.sqrt(2) * np.sqrt(2*degree + 1)
x_min = mean_x - band_width * np.sqrt(cov_x)
x_max = mean_x + band_width * np.sqrt(cov_x)

# Parameters for visualization
n_points_visu, extrema_visu, cov_visu = [1000], [band_width], cov_x
quad_visu = hm.Quad.newton_cotes(n_points_visu, extrema_visu,
                                 mean=[mean_x], cov=[[cov_visu]])
x_visu = quad_visu.discretize('x')
factor_visu = quad_visu.discretize(factor_x)

# }}}
# PLOT HERMITE FUNCTION OF HIGHEST DEGREE {{{

fig, ax = plt.subplots(1, 1)
ax.set_title("Hermite function of degree " + str(degree))
ax.axvline(x=x_min)
ax.axvline(x=x_max)
ax.set_ylim((-2, 2))
h_i = np.zeros(degree + 1)
h_i[degree] = 1
h_i = hm.Series(h_i, mean=[mean_x], cov=[[cov_x]])
# Eh = quad_visu.eval(h_i, degree) * np.sqrt(factor_visu)
Eh = quad_visu.eval(h_i, degree) * factor_visu
ax.plot(x_visu, Eh)
plt.show()

# }}}

# SPECTRAL METHOD FOR STATIONARY EQUATION {{{

mat_operator_x = np.zeros((degree + 1, degree + 1))
for i in range(len(split_operator_x)):
    coeff = split_operator_x[i]
    mat_operator_x += quad_num_x.dvarf(coeff, degree, [0]*i)

mat_operator_y = np.zeros((degree + 1, degree + 1))
for i in range(len(split_operator_y)):
    coeff = split_operator_y[i]
    mat_operator_y += quad_num_y.dvarf(coeff, degree, [0]*i)

tensor_operator_x = hm.tensorize(mat_operator_x, 2, 0)
tensor_operator_y = hm.tensorize(mat_operator_y, 2, 1)
total_op = tensor_operator_x + tensor_operator_y

# Calculate eigenvector in kernel
eigen_values, eigen_vectors = las.eigsh(total_op, k=1, which='SM')
Hu_spec_stat = hm.project(eigen_vectors.T[0], 2, 0)
series_spec_stat = hm.Series(Hu_spec_stat, mean=[mean_x], cov=[[cov_x]])
u_spec_stat_visu = quad_visu.eval(series_spec_stat, degree)

# Comparison between exact solution and solution found using spectral method
fig, ax = plt.subplots(1, 1)
ax.plot(x_visu, abs(u_spec_stat_visu) * factor_visu)
plt.show()

# }}}

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
