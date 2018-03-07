# IMPORT MODULES {{{
import scipy.linalg as la
import sympy as sy
import sympy.printing as syp
import numpy as np
import spectral as sp
import matplotlib.pyplot as plt
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

beta = 1

mean = 0
cov = .1

# potential_p = x**2/2 + 10*sy.cos(x)
potential_p = x**4/4 - x**2/2
potential_q = 0.5*np.log(2*np.pi*cov) + (x - mean)*(x - mean)/(2 * cov)

factor_p = factor_pa.subs(potential_pa, potential_p)
factor_p = factor_p.subs(beta_a, beta)
factor_q = factor_qa.subs(potential_qa, potential_q)
factor_q = factor_q.subs(beta_a, beta)
factor = factor_q*factor_p

multiplication = multiplication_a.subs(potential_pa, potential_p)
multiplication = multiplication.subs(potential_qa, potential_q).doit()
multiplication = multiplication.subs(beta_a, beta).doit()

r_sol = r_sol_a.subs(potential_pa, potential_p)
r_sol = r_sol.subs(beta_a, beta)
v_sol = v_sol_a.subs(potential_pa, potential_p)
v_sol = v_sol.subs(beta_a, beta)
u_sol = u_sol_a.subs(potential_pa, potential_p)
u_sol = u_sol.subs(potential_qa, potential_q)
u_sol = u_sol.subs(beta_a, beta)

# }}}
# DISCRETIZE VARIOUS FUNCTIONS ON GRID {{{

# For numerics
degree = 100
degrees = np.arange(degree + 1)
n_points_num = degree + 1
quad_num = sp.Quad(n_points_num, dim=1, mean=[mean], cov=[[cov]])
x_num = quad_num.discretize('x')
factor_p_num = quad_num.discretize(factor_p)
factor_q_num = quad_num.discretize(factor_q)
factor_num = quad_num.discretize(factor)
multiplication_num = quad_num.discretize(multiplication)
u_sol_num = quad_num.discretize(u_sol)
v_sol_num = quad_num.discretize(v_sol)
r_sol_num = quad_num.discretize(r_sol)

# Calculate limits of resolution
x_min = mean + np.sqrt(2) * np.sqrt(2*degree + 1) * np.sqrt(cov)
x_max = mean - np.sqrt(2) * np.sqrt(2*degree + 1) * np.sqrt(cov)

# For visualization
n_points_visu = 300
cov_visu = cov
quad_visu = sp.Quad(n_points_visu, dim=1, mean=[mean], cov=[[cov_visu]])
x_visu = quad_visu.discretize('x')
factor_p_visu = quad_visu.discretize(factor_p)
factor_q_visu = quad_visu.discretize(factor_q)
factor_visu = quad_visu.discretize(factor)
multiplication_visu = quad_visu.discretize(multiplication)
u_sol_visu = quad_visu.discretize(u_sol)
v_sol_visu = quad_visu.discretize(v_sol)
r_sol_visu = quad_visu.discretize(r_sol)
u_approx_visu = quad_visu.eval(quad_num.transform(u_sol, degree), degree)[0]
v_approx_visu = u_approx_visu * factor_q_visu
r_approx_visu = u_approx_visu * factor_visu


def norm_herm(hermite_coeffs):
    return np.sqrt(np.sum(np.square(hermite_coeffs)))


# Normalization
norm_sol_visu = norm_herm(quad_visu.transform(u_sol_visu, degree))
norm_approx_visu = norm_herm(quad_visu.transform(u_approx_visu, degree))
v_sol_visu = v_sol_visu / norm_sol_visu
r_sol_visu = r_sol_visu / norm_sol_visu
v_approx_visu = v_approx_visu / norm_approx_visu
r_approx_visu = r_approx_visu / norm_approx_visu

Hu_sol = quad_num.transform(u_sol, degree)
Hu_sol = Hu_sol / norm_herm(Hu_sol)

# }}}
# NUMERICAL METHOD

# Plot Hermite function of highest degree
fig, ax = plt.subplots(1, 1)
ax.set_title("Hermite function of degree " + str(degree))
ax.axvline(x=x_min)
ax.axvline(x=x_max)
ax.set_ylim((-2, 2))
h_i = np.zeros(degree + 1)
h_i[degree] = 1
Eh = quad_visu.eval(h_i, degree)[0] * factor_q_visu
ax.plot(x_visu, Eh)
plt.show()

# Time step and number of iterations
dt = 2e-3*cov
n_iter = 100000

# Eigenvalues of the operator
eigenvalues_gaussian = - np.arange(degree + 1)/cov

# Initial condition
u_init = factor_q / factor_p
u_num = quad_num.discretize(u_init)
Hu = quad_num.transform(u_num, degree)

# Exact solution on fine grid
# v_approx_visu = quad_visu.eval(quad_num.transform(u_sol_num, degree), degree)[0]*factor_q_visu
# norm_aux = quad_num.transform(u_exact_num/sqrt_gaussian_num, degree)
# norm = np.sqrt(np.sum(np.square(norm_aux)))
# v_exact_visu = u_exact_visu / norm
# v_approx_visu = u_approx_visu / norm

# plt.plot(x_visu, u_exact_visu)
# plt.show()

# Get matrix representation of operator in Hermite space:
diag_op = np.diag(eigenvalues_gaussian)
multiplication_op = quad_visu.varf(multiplication_visu, degree)
total_op = diag_op + multiplication_op

# eigen_values, eigen_vectors = la.eigs(total_op, k=1, which='SM')

# # Print approximating functions to Schrődinger equation
# fig, ax = plt.subplots(1, 1)
# plt.ion()

# for i in range(10):
#     ax.set_title("First eigenfunctions of Fokker-Planck operator")
#     h_i = np.zeros(degree + 1)
#     h_i[i] = 1
#     Eh = quad_visu.eval(h_i, degree)[0] * factor_q_visu
#     ax.plot(x_visu, Eh)
#     ax.set_ylim((-2, 2))
#     plt.draw()
#     plt.pause(.01)
#     ax.clear()
# plt.close()
# # Plot first eigenvectors



# # Simple integration
# for i in range(n_iter):

#     # Plotting {{{
#     if i % 100 == 0:
#         plt.pause(.01)

#         # Representation of u on fine grid
#         u_visu = quad_visu.eval(Hu, degree)[0]

#         # Plot solution in flat space
#         ax11.clear()
#         ax11.set_title("Solution to Schrődinger equation")
#         ax11.plot(x_visu, u_visu * factor_q_visu)
#         ax11.plot(x_visu, v_sol_visu)
#         ax11.plot(x_visu, v_approx_visu)

#         # Plot solution in real space
#         ax12.clear()
#         ax12.set_title("Solution to Fokker-Planck equation")
#         ax12.plot(x_visu, u_visu * factor_visu)

#         # Plot Hermite transform
#         ax21.clear()
#         ax21.set_title("Hermite coefficients of the solution")
#         ax21.bar(degrees, Hu)

#         plt.draw()
#     # }}}

#     # Normalization
#     Hu = la.expm(total_op)

# Create figure with two subplots
fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)

# Activate interactive plotting
plt.ion()

# Splitting method
for i in range(n_iter):

    # Plotting {{{
    if i % 100 == 0:
        plt.pause(.01)

        # Representation of u on fine grid
        u_visu = quad_visu.eval(Hu, degree)[0]

        # Plot solution in flat space
        ax11.clear()
        ax11.set_title("Solution to Schrődinger equation")
        ax11.plot(x_visu, u_visu * factor_q_visu)
        ax11.plot(x_visu, v_sol_visu)
        ax11.plot(x_visu, v_approx_visu)

        # Plot solution in real space
        ax12.clear()
        ax12.set_title("Solution to Fokker-Planck equation")
        ax12.plot(x_visu, u_visu * factor_visu)

        # Plot Hermite transform
        ax21.clear()
        ax21.set_title("Hermite coefficients of the numerical solution")
        ax21.bar(degrees, Hu)

        # Plot error for Hermite coefficients
        ax22.clear()
        ax22.set_title("Hermite coefficients of the stationary solution")
        ax22.bar(degrees, Hu_sol)

        plt.draw()
    # }}}

    # h_n = quad.transform(u_n, degree)
    Hu = Hu + dt/2 * eigenvalues_gaussian * Hu
    u_num = quad_num.eval(Hu, degree)[0]
    Hu = Hu + dt * quad_num.transform(u_num*multiplication_num, degree)
    Hu = Hu + dt/2 * eigenvalues_gaussian * Hu

    # Normalization
    Hu = Hu/norm_herm(Hu)

    # u_n = u_n + dt * diff_linear_n * u_n

# plt.plot(points[0], potential_n)

import importlib
importlib.reload(sp)
