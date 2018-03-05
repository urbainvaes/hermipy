# IMPORT MODULES {{{
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

potential_p = x**2/2 + 10*sy.cos(x)
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

degree = 100
degrees = np.arange(degree + 1)

n_points_fine = 200
n_points_coarse = degree + 1

quad_fine = sp.Quad(n_points_fine, dim=1, mean=[mean], cov=[[cov]])
quad_coarse = sp.Quad(n_points_coarse, dim=1, mean=[mean], cov=[[cov]])

x_fine = quad_fine.discretize('x')
x_coarse = quad_coarse.discretize('x')

factor_p_fine = quad_fine.discretize(factor_p)
factor_q_fine = quad_fine.discretize(factor_q)
factor_fine = quad_fine.discretize(factor)

factor_p_coarse = quad_coarse.discretize(factor_p)
factor_q_coarse = quad_coarse.discretize(factor_q)
factor_coarse = quad_coarse.discretize(factor)

multiplication_fine = quad_fine.discretize(multiplication)
multiplication_coarse = quad_coarse.discretize(multiplication)

u_sol_fine = quad_fine.discretize(u_sol)
u_sol_coarse = quad_coarse.discretize(u_sol)
v_sol_fine = quad_fine.discretize(v_sol)
v_sol_coarse = quad_coarse.discretize(v_sol)
r_sol_fine = quad_fine.discretize(r_sol)
r_sol_coarse = quad_coarse.discretize(r_sol)

u_approx_fine = quad_fine.eval(quad_coarse.transform(u_sol, degree), degree)[0]
v_approx_fine = u_approx_fine * factor_q_fine
r_approx_fine = u_approx_fine * factor_fine


def norm_herm(hermite_coeffs):
    return np.sqrt(np.sum(np.square(hermite_coeffs)))


# Normalization
norm_sol_fine = norm_herm(quad_fine.transform(u_sol_fine, degree))
norm_approx_fine = norm_herm(quad_fine.transform(u_approx_fine, degree))
v_sol_fine = v_sol_fine / norm_sol_fine
r_sol_fine = r_sol_fine / norm_sol_fine
v_approx_fine = v_approx_fine / norm_approx_fine
r_approx_fine = r_approx_fine / norm_approx_fine

# }}}
# NUMERICAL METHOD

# Print approximating functions to Schrődinger equation
fig, ax = plt.subplots(1, 1)
plt.ion()

for i in [d for d in degrees if d % 2 == 0]:
    ax.set_title("Hermite function of index " + str(i))
    h_i = np.zeros(degree + 1)
    h_i[i] = 1
    Eh = quad_fine.eval(h_i, degree)[0] * factor_q_fine
    ax.plot(x_fine, Eh)
    ax.set_ylim((-2, 2))
    plt.draw()
    plt.pause(.01)
    ax.clear()
plt.close()

# Time step and number of iterations
dt = 2e-3*cov
n_iter = 100000

# Eigenvalues of the operator
eigenvalues_gaussian = - np.arange(degree + 1)/cov

# Initial condition
u_init = factor_q / factor_p
u_coarse = quad_coarse.discretize(u_init)
Hu = quad_coarse.transform(u_coarse, degree)

# Exact solution on fine grid
# v_approx_fine = quad_fine.eval(quad_coarse.transform(u_sol_coarse, degree), degree)[0]*factor_q_fine
# norm_aux = quad_coarse.transform(u_exact_coarse/sqrt_gaussian_coarse, degree)
# norm = np.sqrt(np.sum(np.square(norm_aux)))
# v_exact_fine = u_exact_fine / norm
# v_approx_fine = u_approx_fine / norm

# plt.plot(x_fine, u_exact_fine)
# plt.show()

# Create figure with two subplots
fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)

# Activate interactive plotting
plt.ion()

for i in range(n_iter):

    # Plotting {{{
    if i % 100 == 0:
        plt.pause(.01)

        # Representation of u on fine grid
        u_fine = quad_fine.eval(Hu, degree)[0]

        # Plot solution in flat space
        ax11.clear()
        ax11.set_title("Solution to Schrődinger equation")
        ax11.plot(x_fine, u_fine * factor_q_fine)
        ax11.plot(x_fine, v_sol_fine)
        ax11.plot(x_fine, v_approx_fine)

        # Plot solution in real space
        ax12.clear()
        ax12.set_title("Solution to Fokker-Planck equation")
        ax12.plot(x_fine, u_fine * factor_fine)

        # Plot Hermite transform
        ax21.clear()
        ax21.set_title("Hermite coefficients of the solution")
        ax21.bar(degrees, Hu)

        plt.draw()
    # }}}

    # h_n = quad.transform(u_n, degree)
    Hu = Hu + dt/2 * eigenvalues_gaussian * Hu
    u_coarse = quad_coarse.eval(Hu, degree)[0]
    Hu = Hu + dt * quad_coarse.transform(u_coarse*multiplication_coarse, degree)
    Hu = Hu + dt/2 * eigenvalues_gaussian * Hu

    # Normalization
    Hu = Hu/norm_herm(Hu)

    # u_n = u_n + dt * diff_linear_n * u_n

# plt.plot(points[0], potential_n)

import importlib
importlib.reload(sp)
