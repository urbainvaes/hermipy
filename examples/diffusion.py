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

import os
import math
import argparse
import sympy as sym
import numpy as np
import numpy.linalg as la
import matplotlib
import hermipy as hm
import hermipy.equations as eq
import scipy.integrate as integrate

hm.settings['tensorize'] = True
hm.settings['sparse'] = True
hm.settings['trails'] = False
hm.settings['cache'] = True

# Process arguments
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model', type=str)
parser.add_argument('-dir', '--directory', type=str)
parser.add_argument('-i', '--interactive', action='store_true')
parser.add_argument('-a', '--alpha', type=str)
parser.add_argument('-b', '--beta', type=str)
parser.add_argument('-e', '--epsilon', type=str)
parser.add_argument('-g', '--gamma', type=str)
parser.add_argument('-l', '--lamda', type=str)
parser.add_argument('-s', '--sigma', type=str)
parser.add_argument('-d', '--degree', type=int)
parser.add_argument('-n', '--npoints', type=int)
parser.add_argument('-do', '--overdamped', action='store_true')
parser.add_argument('-d0', '--diff0', action='store_true')
parser.add_argument('-d1', '--diff1', action='store_true')
parser.add_argument('-d2', '--diff2', action='store_true')
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-tg', '--test_gammas', action='store_true')
args = parser.parse_args()

# Directory for output files
dir = ""
if args.directory:
    dir = args.directory + "/"
    os.makedirs(dir, exist_ok=True)

# Matplotlib configuration
matplotlib.rc('font', size=14)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)

if 'DISPLAY' not in os.environ:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

# Equation parameters
r, s = sym.Rational, sym.symbols
e = eq.Generalized_Langevin
x, y, z, w = e.x, e.y, e.z, e.w


def set_param(arg, default, symbol):
    if arg and r(arg) < 0:
        return s(symbol)
    if arg:
        return r(arg)
    return default


ε = set_param(args.epsilon, r(1), 'ε')
β = set_param(args.beta, r(1), 'β')
γ = set_param(args.gamma, r(1), 'γ')
σ = set_param(args.sigma, r(3), 'σ')
Vx = 1 - sym.cos(x)

# Numerical parameters
index_set = 'cube'
s2y, s2z, s2w = r(1, σ), r(1, σ), r(1, σ)
degree = args.degree if args.degree else 10
npoints = args.npoints if args.npoints else 2*degree + 1
kwargs0 = {'degree': degree, 'index_set': index_set}

# Calculate factors {{{
Vqy = sym.Rational(1/2)*y*y/s2y
Vqz = sym.Rational(1/2)*z*z/s2z
Vqw = sym.Rational(1/2)*w*w/s2w

# Normalization constants
zx = integrate.quad(sym.lambdify(x, sym.exp(-β*Vx)), -math.pi, math.pi)[0]
zy = integrate.quad(sym.lambdify(y, sym.exp(-β*y*y/2)), -10, 10)[0]
zz = integrate.quad(sym.lambdify(z, sym.exp(-β*z*z/2)), -10, 10)[0]
zw = integrate.quad(sym.lambdify(w, sym.exp(-β*w*w/2)), -10, 10)[0]

# Map to appropriate space
zfy = sym.sqrt(2*sym.pi*s2y)
zfz = sym.sqrt(2*sym.pi*s2z)
zfw = sym.sqrt(2*sym.pi*s2w)

factor_x = sym.exp(β/2*Vx)  # Fourier
factor_y = sym.exp(- 1/2*Vqy + β/2*y*y/2)  # Hermite functions
factor_z = sym.exp(- 1/2*Vqz + β/2*z*z/2)  # Hermite functions
factor_w = sym.exp(- 1/2*Vqw + β/2*w*w/2)  # Hermite functions
factor_1 = factor_y * factor_z             # 1 extra process
factor_2 = factor_y * factor_z * factor_w  # 2 extra processes
# }}}

# Definition of quadratures
new_q = hm.Quad.gauss_hermite
cov_0, cov_1 = [[s2y]], [[s2y, 0], [0, s2z]]
cov_2 = [[s2y, 0, 0], [0, s2z, 0], [0, 0, s2w]]
args_f = {'dirs': [0], 'factor': factor_x}
args_0 = {'dirs': [1], 'mean': [0]*1, 'cov': cov_0, 'factor': factor_y}
args_1 = {'dirs': [1, 2], 'mean': [0]*2, 'cov': cov_1, 'factor': factor_1}
args_2 = {'dirs': [1, 2, 3], 'mean': [0]*3, 'cov': cov_2, 'factor': factor_2}
qx = hm.Quad.fourier(npoints, **args_f)
q0 = hm.Quad.gauss_hermite(npoints, **args_0)
q1 = hm.Quad.gauss_hermite(npoints, **args_1)
q2 = hm.Quad.gauss_hermite(npoints, **args_2)
# qy_coarse = hm.Quad.fourier(degree + 1, **args_f)
# qxz_coarse = new_q(degree + 1, **args_1)
# qxzw_coarse = new_q(degree + 1, **args_2)
qy_visu = hm.Quad.newton_cotes(n_points=[npoints], extrema=[5],
                               dirs=[1], factor=factor_y)


def diffo():
    eq = hm.equations.Overdamped
    factor = factor_x*sym.sqrt(zx)
    quad = hm.Quad.fourier(npoints, dirs=[0], factor=factor)
    f, backward = eq.f, eq.backward({'β': β, 'V': Vx})
    operator = quad.discretize_op(backward, **kwargs0)
    rhs = quad.transform(-Vx.diff(x), **kwargs0)
    solution = (-operator).solve(rhs)
    Ix = quad.transform(1, **kwargs0)
    diff = quad.discretize_op(f.diff(x), **kwargs0)
    aux = quad.transform(β*Vx.diff(x), **kwargs0)
    # quad.plot(Ix)
    # quad.plot(solution)
    # quad.plot(diff(solution))  # THIS IS NOT CORRECT
    diffusion1 = (1/β) + float(solution*rhs) + (2/β)*float(Ix*diff(solution))
    diffusion2 = (1/β) + float(solution*rhs) + (2/β)*float(aux*solution)
    assert abs(diffusion1 - diffusion2) < 1e-8
    print("Overdamped Langevin: {}".format(diffusion1))
    return diffusion1


# Calculation of the diffusion coefficient with 0 extra process
def diff0():
    quad, params = (qx*q0), {'β': β, 'γ': γ, 'Vx': Vx}
    backward0 = hm.equations.Generalized_Langevin.backward(params)
    operator = quad.discretize_op(backward0, **kwargs0)
    rhs = quad.transform('y', **kwargs0)
    one = quad.transform(1, **kwargs0)

    solution = (- operator).solve(rhs)
    diffusion = float(solution*rhs) * float(zfy/(zx*zy))
    print("With 0 extra process: {}".format(diffusion))

    solution = (- operator).solve(rhs, remove_vec=one)
    diffusion = float(solution*rhs) * float(zfy/(zx*zy))
    print("With 0 extra process: {}".format(diffusion))

    if args.interactive:
        fig, ax = plt.subplots()
        plot = (qx*qy_visu).plot(solution, ax=ax, contours=20)
        plt.colorbar(plot, ax=ax)
        ax.set_xlabel('q')
        ax.set_ylabel('p')
        ax.set_title('$\\gamma = {}$'.format(γ))
        plt.savefig("solution.eps", bbox_inches='tight')
        plt.show()

        # quad.plot(solution, ax=ax, vmin=-200, vmax=200)

    # Calculate autocorrelation function
    # def dudt(t, y):
    #     print("t: {}".format(t))
    #     return_vector = operator.matrix.dot(y)
    #     return return_vector

    # time = np.linspace(0, 10, 101)
    # result = integrate.solve_ivp(dudt, [time[0], time[-1]], rhs.coeffs,
    #                              'RK45', t_eval=time, max_step=.01)
    # vac = np.zeros(len(time))
    # for i, u in enumerate(result.y.T):
    #     u_series = quad.series(u, index_set=index_set)
    #     vac[i] = float(rhs*u_series) * float(zfy/(zx*zy))

    # fig, ax = plt.subplots()
    # ax.plot(time, vac)
    # plt.show()

    return diffusion


# Calculation of the diffusion coefficient with one extra process
def diff1():
    α = set_param(args.alpha, r(1), 'α') / ε**2 / γ
    λ = set_param(args.lamda, r(1), 'λ') / ε
    params, quad = {'β': β, 'α': α, 'λ': λ, 'Vx': Vx}, qx*q1
    backward1 = e.backward(params)
    operator = quad.discretize_op(backward1, **kwargs0)
    rhs = quad.transform('y', **kwargs0)
    solution = (- operator).solve(rhs)
    diffusion = float(solution*rhs) * float(zfy*zfz/(zx*zy*zz))
    print("With 1 extra process: {}".format(diffusion))

    # Calculate autocorrelation function
    # def dudt(t, y):
    #     print("t: {}".format(t))
    #     return_vector = operator.matrix.dot(y)
    #     return return_vector

    # time = np.linspace(0, 30, 201)
    # result = integrate.solve_ivp(dudt, [time[0], time[-1]], rhs.coeffs,
    #                              'RK45', t_eval=time, max_step=.01)
    # vac = np.zeros(len(time))
    # for i, u in enumerate(result.y.T):
    #     u_series = quad.series(u, index_set=index_set)
    #     vac[i] = float(rhs*u_series) * float(zfy*zfz/(zx*zy*zz))

    # fig, ax = plt.subplots()
    # ax.plot(time, vac)
    # plt.show()

    return diffusion


# Calculation of the diffusion coefficient with two extra processes
def diff2(x0=None):
    λ1 = set_param(args.lamda, r(1), 'λ') / ε
    A12 = -1/ε**2/γ
    A21 = 1/ε**2/γ
    A22 = set_param(args.alpha, r(1), 'α') / ε**2 / γ
    params = {'β': β, 'Vx': Vx, 'λ1': λ1, 'λ2': 0,
              'A11': 0, 'A12': A12, 'A21': A21, 'A22': A22}
    backward = e.backward(params)
    qf = qx * q2
    rhs = qf.transform('y', **kwargs0)
    operator = qf.discretize_op(backward, **kwargs0)
    solution = (-operator).solve(rhs, use_gmres=True, x0=x0, tol=5e-2,
                                 callback=lambda rk: print(la.norm(rk)))
    diffusion = float(solution*rhs) * float(zfy*zfz*zfw/(zx*zy*zz*zw))
    print("With 2 extra processes: {}".format(diffusion))
    return diffusion, solution.coeffs


γs = np.logspace(-3, 3, 100) if args.test_gammas else [γ]
diffusion_coefficients = []
do = diffo() if args.overdamped else None
guess = None
for γ in np.flip(γs):
    print("Value of γ: {}".format(γ))
    d0 = diff0() if args.diff0 else None
    d1 = diff1() if args.diff1 else None
    d2, guess = diff2(guess) if args.diff2 else (None, None)
    diffusion_coefficients = []
exit(0)

# # Projections
# qx, qy, qz = quad.project(0), quad.project(1), quad.project(2)
# wx = qx.factor * qx.factor / qx.position.weight()
# wy = qy.factor * qy.factor / qy.position.weight()
# wz = qz.factor * qz.factor / qz.position.weight()
# qxy, qxz, qyz = quad.project([0, 1]), quad.project([0, 2]), quad.project([1, 2])

# # Integral operators
# Ix = qx.transform(wx, degree=degree, index_set=index_set)
# Iy = qy.transform(wy, degree=degree, index_set=index_set)
# Iz = qz.transform(wz, degree=degree, index_set=index_set)
