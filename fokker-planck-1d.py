# IMPORT MODULES {{{
import numpy.linalg as la
import scipy.sparse.linalg as las
import sympy as sy
import sympy.printing as syp
import numpy as np
import matplotlib.pyplot as plt

from libhermite import hermite_python as sp
sy.init_printing()
# }}}
# ABSTRACT SYMBOLIC CALCULATIONS {{{

# Number of space dimensions
dim = 1

# Space variable, time variable and inverse temperature
x, t, beta_a = sy.symbols('x t beta')

# Potential of interest
potential_pa = sy.Function('V')(x)

# Quadratic potential used for the approximation
potential_qa = sy.Function('Vq')(x)

# Solution of the Fokker-Planck equation
r = sy.Function('ρ')(x, t)

# Mapped solution: u = ρ * e^{V/2} * e^{Vq/2} to equation with BK operator
u = sy.Function('u')(x, t)
u_x = sy.Function('u')(x)

# Mapped solution of Schrodinger equation
v = sy.Function('v')(x, t)


# Backward Kolmogorov operator associated with potential
def backward(potential, f):
    dx_potential = sy.diff(potential, x)
    dx_f = sy.diff(f, x)
    dxx_f = sy.diff(dx_f, x)
    return - dx_potential * dx_f + (1/beta_a) * dxx_f


# Fokker-Planck operator associated with potential
def forward(potential, f):
    dx_potential = sy.diff(potential, x)
    dxx_potential = sy.diff(dx_potential, x)
    dx_f = sy.diff(f, x)
    dxx_f = sy.diff(dx_f, x)
    return (dx_potential * dx_f + dxx_potential * f) + (1/beta_a) * dxx_f


# Factors to map Fokker-Planck to backward Kolmogorov operator
factor_pa = sy.exp(-beta_a * potential_pa / 2)
factor_qa = sy.exp(-beta_a * potential_qa / 2)

# Mapping to L² space weighted by Gaussian
factor_a = factor_pa * factor_qa

# Solutions
r_sol_a = sy.exp(-beta_a * potential_pa)
v_sol_a = r_sol_a/factor_pa
u_sol_a = r_sol_a/factor_a

# Fokker-Planck equation to solve (= 0)
fk = sy.diff(r, t) - forward(potential_pa, r)

# Schrodinger equation associated with fk
sch = sy.simplify(fk.subs(r, v*factor_pa).doit()/factor_pa)

# Backward Kolmogorov-like equation
bk = sy.simplify(fk.subs(r, u*factor_a).doit()/factor_a)

# Obtain RHS (if the left-hand side is the time derivative)
operator_rhs = - bk.subs(u, sy.Function('u')(x)).doit()

# Obtain multiplication operator
multiplication_a = operator_rhs.subs(u_x, 1).doit()

# Remainder of the operator
generator_ou = sy.simplify(operator_rhs - multiplication_a*u_x)

# Assert that `generator_ou` is the generator of an OU process
assert generator_ou == backward(potential_qa, u_x)

# Assert that the solutions found above are correct
assert fk.subs(r, r_sol_a).doit().simplify() == 0
assert sch.subs(v, v_sol_a).doit().simplify() == 0
assert bk.subs(u, u_sol_a).doit().simplify() == 0

# }}}
# PRINT TO STDOUT {{{
print("Fokker-Planck equation to solve: ")
syp.pprint(fk)

print("Mapping to a Schrődinger equation")
syp.pprint(sch)

print("Mapping to an equation with BK operator")
syp.pprint(bk)

print("Operator in the right-hand side")
syp.pprint(operator_rhs)

print("Multiplication part of the operator")
syp.pprint(multiplication_a)

print("Remainder of the operator")
syp.pprint(generator_ou)
# }}}
# EVALUATE ABSTRACT EXPRESSIONS FOR PROBLEM AT HAND {{{

# Parameter in diffusion
beta = 2.

# For numerical approximation
mean = .2
cov = .1

# potential_p = x**2/2 + 10*sy.cos(x)
potential_p = x**4/4 - x**2/2
potential_q = 0.5*np.log(2*np.pi*cov)/beta + (x - mean)*(x - mean)/(2*beta*cov)

factor_p = factor_pa.subs(potential_pa, potential_p)
factor_p = factor_p.subs(beta_a, beta)
factor_q = factor_qa.subs(potential_qa, potential_q)
factor_q = factor_q.subs(beta_a, beta)
factor = factor_q*factor_p

multiplication = multiplication_a.subs(potential_pa, potential_p)
multiplication = multiplication.subs(potential_qa, potential_q).doit()
multiplication = multiplication.subs(beta_a, beta).doit()

u_sol = u_sol_a.subs(potential_pa, potential_p)
u_sol = u_sol.subs(potential_qa, potential_q)
u_sol = u_sol.subs(beta_a, beta)

# Normalization using precise quadrature
n_precise = 200
quad_precise = sp.Quad.gauss_hermite(n_precise, dim=1,
                                     mean=[mean], cov=[[cov]])
norm_sol = np.sqrt(quad_precise.integrate(u_sol*u_sol))
u_sol = u_sol / norm_sol

# Hermite transform of the exact solution
Hu_sol = quad_precise.transform(u_sol, n_precise - 1)

# }}}
# DISCRETIZE VARIOUS FUNCTIONS ON GRID {{{

# For numerics
degree = 50
degrees = np.arange(degree + 1)
n_points_num = degree + 1
quad_num = sp.Quad.gauss_hermite(n_points_num, dim=1, mean=[mean], cov=[[cov]])
multiplication_num = quad_num.discretize(multiplication)

# Calculate limits of resolution
x_min = mean + np.sqrt(2) * np.sqrt(2*degree + 1) * np.sqrt(cov)
x_max = mean - np.sqrt(2) * np.sqrt(2*degree + 1) * np.sqrt(cov)

# Parameters for visualization
n_points_visu = [1000]
extrema_visu = [np.sqrt(2) * np.sqrt(2*degree + 1)]
cov_visu = 2*cov
quad_visu = sp.Quad.newton_cotes(n_points_visu, extrema_visu,
                                 mean=[mean], cov=[[cov_visu]])
x_visu = quad_visu.discretize('x')
u_sol_visu = quad_visu.discretize(u_sol)
factor_q_visu = quad_visu.discretize(factor_q)
factor_visu = quad_visu.discretize(factor)

# }}}
# PLOT HERMITE FUNCTION OF HIGHEST DEGREE {{{

fig, ax = plt.subplots(1, 1)
ax.set_title("Hermite function of degree " + str(degree))
ax.axvline(x=x_min)
ax.axvline(x=x_max)
ax.set_ylim((-2, 2))
h_i = np.zeros(degree + 1)
h_i[degree] = 1
h_i = sp.Series(h_i, mean=[mean], cov=[[cov]])
Eh = quad_visu.eval(h_i, degree) * factor_q_visu
ax.plot(x_visu, Eh)
plt.show()

# }}}
# SPECTRAL METHOD FOR STATIONARY EQUATION {{{

# Eigenvalues of the operator
eigenvalues_gaussian = - np.arange(degree + 1)/cov/beta

# Get matrix representation of operator in Hermite space:
diag_op = np.diag(eigenvalues_gaussian)
multiplication_op = quad_precise.varf(multiplication, degree)
total_op = diag_op + multiplication_op

# Calculate eigenvector in kernel
eigen_values, eigen_vectors = las.eigsh(total_op, k=1, which='SM')
Hu_spec_stat = eigen_vectors.T[0]
series_spec_stat = sp.Series(Hu_spec_stat, mean=[mean], cov=[[cov]])
u_spec_stat_visu = quad_visu.eval(series_spec_stat, degree)

# Comparison between exact solution and solution found using spectral method
fig, ax = plt.subplots(1, 1)
ax.plot(x_visu, abs(u_sol_visu) * factor_q_visu)
ax.plot(x_visu, abs(u_spec_stat_visu) * factor_q_visu)
plt.show()

# }}}

# Time step and number of iterations
dt = 2e-4*cov
n_iter = 10**6

# Initial condition
u_init = factor_q / factor_p
u_num = quad_num.discretize(u_init)
Hu = quad_num.transform(u_num, degree)

# Create figure with two subplots
fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)

# Activate interactive plotting
plt.ion()

# Splitting method
for i in range(n_iter):

    # Plotting {{{
    if i % 1000 == 0:
        plt.pause(.01)

        # Representation of u on fine grid
        u_visu = quad_visu.eval(Hu, degree)

        # Plot solution in flat space
        ax11.clear()
        ax11.set_title("Solution to Schrődinger equation")
        ax11.plot(x_visu, u_visu * factor_q_visu)
        ax11.plot(x_visu, u_sol_visu * factor_q_visu)

        # Plot solution in real space
        ax12.clear()
        ax12.set_title("Solution to Fokker-Planck equation")
        ax12.plot(x_visu, u_visu * factor_visu)

        # Plot Hermite transform
        ax21.clear()
        ax21.set_title("Hermite coefficients of the numerical solution")
        ax21.bar(degrees, Hu.coeffs)

        # Plot error for Hermite coefficients
        ax22.clear()
        ax22.set_title("Hermite coefficients of the stationary solution")
        ax22.bar(degrees, Hu_sol.coeffs[0:degree + 1])

        plt.draw()
    # }}}

    Hu.coeffs = Hu.coeffs + dt/2 * eigenvalues_gaussian * Hu.coeffs
    u_num = quad_num.eval(Hu, degree)
    Hu.coeffs = Hu.coeffs + dt * quad_num.transform(u_num*multiplication_num, degree).coeffs
    Hu.coeffs = Hu.coeffs + dt/2 * eigenvalues_gaussian * Hu.coeffs

    # Normalization
    Hu.coeffs = Hu.coeffs/la.norm(Hu.coeffs, 2)

import importlib
importlib.reload(sp)
