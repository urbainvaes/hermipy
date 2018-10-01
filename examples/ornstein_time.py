#!/usr/bin/env python
#
# Copyright (C) 2018 Urbain Vaes
#
# This file is part of hermipy, a python/C++ library for automating the
# Hermite Galerkin method.
#
# hermipy is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hermipy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# Import modules {{{
import argparse
import scipy.integrate
import sympy as sym
import numpy as np
import hermipy as hm
import hermipy.equations as eq
import matplotlib.pyplot as plt
import matplotlib

sym.init_printing()

# }}}
# Parse options {{{
parser = argparse.ArgumentParser()
parser.add_argument('-tp', '--test_plots', action='store_true')
parser.add_argument('-tc', '--test_convergence', action='store_true')
parser.add_argument('-i', '--interactive', action='store_true')
parser.add_argument('-dir', '--directory', type=str)
parser.add_argument('-m0', '--mass0', type=float)
parser.add_argument('-e', '--epsilon', type=str)
parser.add_argument('-b', '--beta', type=float)
parser.add_argument('-t', '--theta', type=str)
parser.add_argument('-m', '--mass', type=float)
parser.add_argument('-d', '--degree', type=int)
parser.add_argument('--method', type=str)
args = parser.parse_args()

# Set hermipy options
hm.settings['tensorize'] = True
hm.settings['sparse'] = True
hm.settings['trails'] = False
hm.settings['cache'] = True

# Set matplotlib options
matplotlib.rc('font', size=14)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)
# }}}
# Data and parameters for numerical simulation {{{


def set_param(arg, default, symbol):
    r, s = sym.Rational, sym.symbols
    if arg and r(arg) < 0:
        return s(symbol)
    if arg:
        return r(arg)
    return default


# Variables and function
x, y = eq.McKean_Vlasov.x, eq.McKean_Vlasov.y
f = eq.McKean_Vlasov.f

r, s = sym.Rational, sym.symbols
β = set_param(args.beta, r(1), 'β')
ε = set_param(args.epsilon, r(1, 2), 'ε')
θ = set_param(args.theta, 1, 'θ')
m = set_param(args.mass, s('m'), 'm')
m0x, m0y = set_param(args.mass0, 1, 'm0'), 1
Vp, Vy = x*x/2, y*y/2

params = {'β': β, 'ε': ε, 'γ': 0, 'θ': θ, 'm': m, 'Vp': Vp, 'Vy': Vy}
s2x, s2y = r(1, 5), r(1, 5)
Vqx, Vqy = r(1/2)*x*x/s2x, r(1/2)*y*y/s2y
forward = eq.McKean_Vlasov.equation(params)

factor_x = sym.exp(- 1/2 * (Vqx + β*Vp))
factor_y = sym.exp(- 1/2 * (Vqy + Vy))
# factor_x = sym.exp(- 1/2 * Vqx)
# factor_y = sym.exp(- 1/2 * Vqy)
# factor_x = sym.exp(- β*Vp)
# factor_y = sym.exp(- Vy)
factor = factor_x * factor_y
degree, index_set = args.degree if args.degree else 60, 'cube'
kwargs0 = {'degree': degree, 'index_set': index_set}
n_points_num = 2*degree + 1

# Calculation of the solution
new_gh, new_nc = hm.Quad.gauss_hermite, hm.Quad.newton_cotes
cov = [[s2x, 0], [0, s2y]]
quad = new_gh(n_points_num, dim=2, mean=[0]*2, cov=cov, factor=factor)
quad_visu = new_nc([200, 200], [4, 4], factor=factor)

# Operators
m_operator = forward.diff(m)
r_operator = (forward - m*m_operator).cancel()
m_mat = quad.discretize_op(m_operator, **kwargs0)
r_mat = quad.discretize_op(r_operator, **kwargs0)

# Calculate projections
qx, qy = quad.project(0), quad.project(1)
wx = qx.factor * qx.factor / qx.position.weight()
wy = qy.factor * qy.factor / qy.position.weight()

# Integral and moment operators
Ix = qx.transform(wx, degree=degree, index_set=index_set)
Iy = qy.transform(wy, degree=degree, index_set=index_set)
mx1 = qx.transform(wx * x, degree=degree, index_set=index_set)
my1 = qy.transform(wy * y, degree=degree, index_set=index_set)

# Time
time = np.linspace(0, 3, 301)
# }}}
# Exact solution {{{
dimension = 2
μ0, Σ0 = sym.Matrix([m0x, m0y]), sym.Matrix([[1, 0], [0, 1]])
B = sym.Matrix([[-1, 1/sym.sqrt(β)/ε], [0, -1/ε**2]])
K = sym.Matrix([[- θ, 0], [0, 0]])
D = sym.Matrix([[0, 0], [0, 1/ε**2]])
s, t = sym.symbols('s t')
μt = sym.exp(B*t)*μ0
Σt = sym.exp(t*(B + K))*Σ0*sym.exp(t*(B + K).T) + \
    2*(sym.exp(s*(B + K))*D*sym.exp(s*(B + K).T)).integrate((s, 0, t))
coords = sym.Matrix([x, y]) - μt
ρt = 1/sym.sqrt((2*sym.pi)**dimension*Σt.det()) \
        * sym.exp(-1/2*(coords.T*Σt.inv()*coords)[0, 0])
# }}}
# Solution of the problem {{{


def solve(method, subdegree=degree):

    # Subdegree
    sub_r_mat = r_mat.subdegree(subdegree)
    sub_m_mat = m_mat.subdegree(subdegree)
    sub_Ix = Ix.subdegree(subdegree)
    sub_Iy = Iy.subdegree(subdegree)
    sub_mx1 = mx1.subdegree(subdegree)

    # Initial condition
    u = sym.exp(-(x-m0x)*(x-m0x)/2) * sym.exp(-(y - m0y)*(y - m0y)/2)
    u = u / quad.integrate(u, flat=True)
    t = quad.transform(u, degree=subdegree, index_set=index_set)

    if method == "ode45":

        def dfdt(t, y):
            series = quad.series(y, index_set=index_set)
            m = float(sub_mx1*(sub_Iy*series)) / float(sub_Ix*(sub_Iy*series))
            print("t: {}, m: {}".format(t, m))
            mat = sub_r_mat.matrix + m*sub_m_mat.matrix
            return_vector = mat.dot(y)
            return return_vector

        result = scipy.integrate.solve_ivp(dfdt, [time[0], time[-1]], t.coeffs,
                                           'RK45', t_eval=time, max_step=.01)
        result = [quad.series(y, index_set=index_set) for y in result.y.T]

    elif method == "semi_explicit":

        m = float(sub_mx1*(sub_Iy*t))
        eye = quad.varf('1', degree=subdegree, index_set=index_set)
        steps = np.diff(time)
        scheme = 'crank'
        result = [t]

        for i, dt in enumerate(steps):

            m = float(sub_mx1*(sub_Iy*t)) / float(sub_Ix*(sub_Iy*t))
            print("i: {}, t: {}, m: {}".format(i, time[i], m))

            mat = sub_r_mat + m*sub_m_mat
            if scheme == 'backward':
                total_op = eye - dt*mat
                new_t = total_op.solve(t)

            if scheme == 'crank':
                crank_left = eye - dt*mat*(1/2)
                crank_right = eye + dt*mat*(1/2)
                new_t = crank_left.solve(crank_right(t))

            # t = new_t / float(In*new_t)
            t = new_t
            result.append(t)

    return result


# }}}
# Convergence {{{
def error_Linf(result, prefix=""):
    errors = np.zeros(len(result))
    for i, t in enumerate(result):
        exact = (ρt.subs(sym.symbols('t'), time[i]))
        t_exact = quad.transform(exact, degree=degree, index_set=index_set)
        errors[i] = quad.norm(t.subdegree(degree) - t_exact, n=1, flat=True)
        print('{}: i: {}, Δ: {}'.format(prefix, i, errors[i]))
    return np.max(errors)


if args.test_convergence:
    degrees = list(range(20, degree, 5))
    errors = np.zeros(len(degrees))
    method = args.method if args.method else "ode45"

    for i, d in enumerate(degrees):
        result = solve(method, subdegree=d)
        errors[i] = error_Linf(result)
    np.save("comp_exact_gaussian_time_degrees_{}".format(method), degrees)
    np.save("comp_exact_gaussian_time_errors_{}".format(method), errors)

# }}}
# Plots {{{
if args.test_plots:

    result = solve('semi_explicit')

    plt.ion()
    fig, ax = plt.subplots(2)

    def plot(i, t, title=None):

        ax[0].clear()
        quad_visu.plot(t, ax=ax[0], bounds=False, title=title)
        # ax[0].set_xlim(-4, 4)

        ax[1].clear()
        exact = quad_visu.discretize(ρt.subs(sym.symbols('t'), time[i]))
        quad_visu.plot(exact - quad_visu.eval(t), ax=ax[1], bounds=False)

        plt.draw()
        plt.pause(0.01)

    for i, t in enumerate(result):
        if i % 100 == 0:
            plot(i, t)
        exact = (ρt.subs(sym.symbols('t'), time[i]))
        t_exact = quad.transform(exact, degree=degree, index_set=index_set)
        Δ = quad.norm(t - t_exact, n=1, flat=True)
        print('i: {}, Δ: {}'.format(i, Δ))

    plot_times = [0, 1, 2, 3]
    plt.ioff()
    for t in plot_times:
        i = np.argmin(np.abs(np.asarray(time) - t))
        fig, ax = plt.subplots()
        quad_visu.plot(result[i], ax=ax, bounds=False,
                       title="Time: " + str(time[i]))
        plt.savefig('solution_time_dependent_gaussian_case_time={}.pdf'.format(t),
                    bbox_inches='tight')
        plt.show()

# }}}
