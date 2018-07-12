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


#  TODO: Error with project (urbain, Fri 27 Apr 2018 05:23:07 PM BST)

# IMPORT MODULES {{{
import multiprocessing

import argparse
import sympy as sym
# import sympy.plotting as splot
import importlib
import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as las

import hermipy
import hermipy.stats
import hermipy.equations as eq
import hermipy.core as core

from sympy.printing import pprint

sym.init_printing()

# }}}
# PARSE OPTIONS {{{
parser = argparse.ArgumentParser()

parser.add_argument('-tc', '--convergence', action='store_true',
                    help='Run convergence test')
parser.add_argument('-tb', '--bifurcation', action='store_true',
                    help='Run bifurcation test')

parser.add_argument('-c', '--config', type=str, help='Configuration file')
parser.add_argument('-b', '--beta', type=str, help='Value of β')
parser.add_argument('-t', '--theta', type=str, help='Value of θ')
parser.add_argument('-e', '--epsilon', type=str, help='Value of ε')
parser.add_argument('-g', '--gamma', type=str, help='Value of γ')
parser.add_argument('-m', '--mass', type=str, help='Value of m')

parser.add_argument('-p', '--plots', action='store_true',
                    help='Enable plots')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Enable verbose output')
parser.add_argument('--cache', action='store_true',
                    help='Enable verbose output')
args = parser.parse_args()


if args.config:
    config = importlib.import_module(args.config)
else:
    config = importlib.import_module("examples.config")

if args.beta:
    config.eq['β'] = sym.Rational(args.beta)
if args.theta:
    config.eq['θ'] = sym.Rational(args.theta)
if args.epsilon:
    config.eq['ε'] = sym.Rational(args.epsilon)
if args.gamma:
    config.eq['γ'] = sym.Rational(args.gamma)
if args.verbose:
    config.misc['verbose'] = args.verbose

if args.plots:
    config.misc['plots'] = args.plots
elif 'plots' not in config.misc:
    config.misc['plots'] = False

if args.cache:
    config.misc['cache'] = args.cache
elif 'plots' not in config.misc:
    config.misc['cache'] = False

if args.verbose:
    config.misc['verbose'] = args.verbose
elif 'verbose' not in config.misc:
    config.misc['verbose'] = False


def vprint(*args, **kwargs):
    if config.misc['verbose']:
        print(*args, **kwargs)


# Set library option
hermipy.settings.update(config.misc)

if config.misc['plots']:
    import matplotlib.pyplot as plt
    # import examples.plot as plots

# }}}
# DATA AND PARAMETERS FOR NUMERICAL SIMULATION {{{
vprint("Symbolic calculations")

# Variables and function
x, y = eq.McKean_Vlasov.x, eq.McKean_Vlasov.y
f = eq.McKean_Vlasov.f

# Parameters
params = {}


class Parameter():

    def __init__(self, symbol, value, params={}, type="main"):
        self.symbol = symbol
        self.value = value
        self.params = params
        self.type = type

    def eval(self):
        value = self.value
        sym_val = {params[k].symbol: params[k].value for k in params}
        m, β = params['m'].symbol, params['β'].symbol
        while value.free_symbols - {x, y, m, β} != set():
            for symbol in self.value.free_symbols:
                if symbol in sym_val:
                    value = value.subs(symbol, sym_val[symbol])
        # if value.free_symbols == set():
        #     value = float(value)
        return value


def eq_params():
    eq_params = {}
    for key, symbol in eq.McKean_Vlasov.params().items():
        value = config.eq[key] if key in config.eq else symbol
        eq_params[key] = Parameter(symbol, value, type='equation')
    extra_symbols = config.eq['Vp'].free_symbols - {x}
    for symbol in extra_symbols:
        key = str(symbol)
        eq_params[key] = Parameter(symbol, config.eq[key], type='equation')
    return eq_params


def vy():
    Vy = Parameter(sym.Function('Vy')(y), y*y/2, type='equation')
    return {'Vy': Vy}


def vq():
    r = sym.Rational
    symb = sym.symbols
    β = params['β']
    μx = Parameter(symb('μx', real=True), config.num['μx'])
    σx = Parameter(symb('σx', real=True, positive=True), config.num['σx'])
    Vq = Parameter(sym.Function('Vq')(x), r(1, 2) *
                   (x - μx.symbol)*(x - μx.symbol)/(σx.symbol*β.symbol))
    return {'μx': μx, 'σx': σx, 'Vq': Vq}


# Forward operator
def forward_op(symbolic):
    eq_params = {}
    for k in params:
        if params[k].type != 'equation':
            continue
        if symbolic == 2:
            value = params[k].symbol
        elif symbolic == 1:
            is_func = params[k].value.free_symbols != set()
            value = params[k].value if is_func else params[k].symbol
        elif symbolic == 0:
            value = params[k].eval()
        eq_params[k] = value
    return eq.McKean_Vlasov.equation(eq_params)


def factors(symbolic, λ):
    β = params['β'].value if symbolic == 0 else params['β'].symbol
    if symbolic == 2:
        Vy = params['Vy'].symbol
        Vp = params['Vp'].symbol
        Vq = params['Vq'].symbol
    elif symbolic == 1:
        Vy = params['Vy'].value
        Vp = params['Vp'].value
        Vq = params['Vq'].value
    elif symbolic == 0:
        Vy = params['Vy'].eval()
        Vp = params['Vp'].eval()
        Vq = params['Vq'].eval()
    factor_x = sym.exp(-β*(λ*Vq + (1-λ)*Vp))
    factor_y = sym.exp(-Vy)
    factor = factor_x * factor_y
    return factor_x, factor_y, factor


def evaluate(expressions):
    subs_list = [(params[key].symbol, params[key].eval())
                 for key in params]
    result = []
    for expr in expressions:
        result.append(expr.subs(subs_list).doit().expand())
    return result


# Parameters
params.update(eq_params())
params.update({**vy(), **vq()})

# Forward operator
forward = forward_op(config.misc['symbolic'])

# Mapped backward operator
factor_x, factor_y, factor = factors(config.misc['symbolic'], config.num['λ'])
backward = eq.map_operator(forward, f, factor)

forward, backward = evaluate([forward, backward])
factor_x, factor_y, factor = evaluate([factor_x, factor_y, factor])

if config.misc['verbose']:
    pprint(forward)
    pprint(backward)

for key in params:
    params[key].value = params[key].eval()

degree = config.num['degree']
n_points_num = config.num['n_points_num']

# }}}
# DEFINE QUADRATURES {{{
vprint("Defining quadratures")


def compute_quads():
    μx = params['μx'].value
    σx = params['σx'].value
    σy = 1

    quad_num = hermipy.Quad.gauss_hermite(
            config.num['n_points_num'], dim=2,
            mean=[μx, 0], cov=[[σx, 0], [0, σy]])

    # Calculate limits of resolution
    # band_width = np.sqrt(2) * np.sqrt(2*degree + 1)
    # x_min = μx - band_width * np.sqrt(σx)
    # x_max = μx + band_width * np.sqrt(σx)
    # y_min = 0 - band_width * np.sqrt(σy)
    # y_max = 0 + band_width * np.sqrt(σy)

    # For visualization
    nv = 200
    # bounds_x, bounds_y = band_width*np.sqrt('σx'), band_width*np.sqrt(cov_y)
    # bounds_x, bounds_y = 5*np.sqrt('σx'), 5*np.sqrt(cov_y)
    bounds_x, bounds_y = 3, 3
    quad_visu = hermipy.Quad.newton_cotes([nv, nv], [bounds_x, bounds_y])
    return quad_num, quad_visu


quad_num, quad_visu = compute_quads()

# }}}
# SPECTRAL METHOD FOR STATIONARY EQUATION {{{

vprint("Splitting operator in m-linear and m-independent parts")
m_operator = backward.diff(params['m'].symbol)
r_operator = (backward - params['m'].symbol*m_operator).cancel()

vprint("Discretizing operators")
m_mat = quad_num.discretize_op(m_operator, f, degree, 2)
r_mat = quad_num.discretize_op(r_operator, f, degree, 2)

vprint("Solving eigenvalue problem")
eig_vals, eig_vecs = r_mat.eigs(k=4, which='LR')
eig_vals = list(reversed(eig_vals))
eig_vecs.reverse()

# if config.misc['verbose']:
#     hermite.stats.print_stats()

# if args.mass:
#     print(compute_moment1(float(args.mass)))
#     matrix = r_mat + float(args.mass) * m_mat
#     eig_vals, eig_vecs = cache(las.eigs)(matrix, k=1, which='LR')
#     ground_state = np.real(eig_vecs.T[0])
#     ground_state = ground_state * np.sign(ground_state[0])
#     ground_series = quad_num.series(ground_state)
#     ground_state_eval = quad_num.eval(ground_series)
#     factor_eval = quad_num.discretize(factor)
#     ground_state_eval = ground_state_eval * factor_eval
#     norm = quad_num.integrate(ground_state_eval, l2=True)
#     ground_state_eval = ground_state_eval / norm
#     fig, ax = plt.subplots(1, 1)
#     cont = quad_visu.plot(ground_series, degree, factor, ax)
#     plt.colorbar(cont, ax=ax)
#     plt.show()

def plot():

    print("[plot]")

    def plot_eigenfunctions():
        fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
        axes = (ax11, ax12, ax21, ax22)
        for e_val, series, ax in zip(eig_vals, eig_vecs, axes):
            series = series * np.sign(series.coeffs[0])
            cont = quad_visu.plot(series, factor, ax=ax, bounds=False)
            ax.set_title("Eigenvalue: " + str(e_val))
            plt.colorbar(cont, ax=ax)
        plt.show()

    def plot_hermite_functions():

        def plot_hf(degrees, quad_num, quad_visu, ax):
            dim = len(degrees)
            adim_width = [np.sqrt(2) * np.sqrt(2*d + 1) for d in degrees]
            pos, bounds = quad_num.position, []
            for i in range(dim):
                width = adim_width[i] * np.sqrt(pos.cov[i][i])
                bounds.append([pos.mean[i] - width, pos.mean[i] + width])
            if dim >= 1:
                ax.axvline(x=bounds[0][0])
                ax.axvline(x=bounds[0][1])
            if dim == 2:
                ax.axhline(y=bounds[1][0])
                ax.axhline(y=bounds[1][1])
            deg_max = sum(degrees)
            index_set = "triangle"
            hf = np.zeros(core.iterator_size(2, deg_max, index_set=index_set))
            hf[core.iterator_index(degrees, index_set="triangle")] = 1
            hf = quad_num.series(hf, quad_num.position)
            factor = quad_num.position.weight()
            return quad_visu.plot(hf, factor, ax, bounds=False)

        fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
        plot_hf((degree, 0), quad_num, quad_visu, ax11)
        plot_hf((0, degree), quad_num, quad_visu, ax12)
        plot_hf((degree // 2, degree // 2), quad_num, quad_visu, ax21)
        plt.show()

    def plot_ground_state():
        ground_state = eig_vecs[0] * np.sign(eig_vecs[0].coeffs[0])
        factors = {'x': factor_x, 'y': factor_y}

        # Plot the projections
        degrees = np.arange(degree + 1)
        series_x = ground_state.project(0)
        series_y = ground_state.project(1)
        quad_x = quad_visu.project(0)
        quad_y = quad_visu.project(1)
        fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
        quad_x.plot(series_x, factors['x'], ax11)
        ax11.set_title("Marginal distribution in x direction")
        ax12.bar(degrees, series_x.coeffs)
        ax12.set_title("Hermite coefficients of x marginal")
        quad_y.plot(series_y, factors['y'], ax21)
        ax21.set_title("Marginal distribution in y direction")
        ax22.bar(degrees, series_y.coeffs)
        ax22.set_title("Hermite coefficients of y marginal")
        plt.show()

    def plot_comparison_with_asym():
        fig, ax = plt.subplots(1, 1)
        as_sol = sym.exp(- params['β'].value * params['Vp'].eval())/factor_x
        as_series = quad_num.project(0).transform(as_sol, degree, norm=True)
        quad_visu.project(0).plot(as_series, factor_x, ax)
        plt.show()

    def plot_discretization_error():
        fig, ax = plt.subplots(1, 2)
        quad_visu_x = quad_visu.project(0)
        as_sol = sym.exp(- params['β'].value * params['Vp'].eval())
        as_series = quad_num.project(0).transform(as_sol/factor_x, degree)
        norm_as = la.norm(as_series.coeffs, 2)
        as_series.coeffs = as_series.coeffs / norm_as
        quad_visu_x.plot(as_series, degree, factor_x, ax[0])
        sol_x = quad_visu_x.discretize(as_sol)
        ax[0].plot(quad_visu_x.discretize('x'), sol_x/norm_as)
        plt.show()

    plot_eigenfunctions()
    plot_hermite_functions()
    plot_ground_state()
    plot_comparison_with_asym()


plot()


def bifurcation(factor_x, m_init):

    def compute_m(series):
        quad_x = quad_num.project(0)
        series_x = series.project(0)
        x_d = quad_x.discretize('x')
        factor_d = quad_x.discretize(factor_x)
        solution_d = quad_x.eval(series_x) * factor_d
        weight_d = quad_x.discretize(quad_x.position.weight())
        moment0 = quad_x.integrate(solution_d/weight_d)
        moment1 = quad_x.integrate(solution_d/weight_d*x_d)
        return moment1 / moment0

    def compute_moment1(m_val):
        if config.misc['verbose']:
            print("Solving the eigenvalue problem for m = " + str(m_val) + ".")
        mat_operator = r_mat + m_mat * m_val
        eig_vals, eig_vecs = mat_operator.eigs(k=1, which='LR')
        ground_state = eig_vecs[0] * np.sign(eig_vecs[0].coeffs[0])
        value = compute_m(ground_state)
        if config.misc['verbose']:
            print(value)
        if config.misc['plots']:
            fig, ax = plt.subplots(1, 1)
            cont = quad_visu.plot(ground_state, factor, ax)
            plt.colorbar(cont, ax=ax)
            plt.show()
        return value

    # m_values = np.linspace(-2, 2, 41)
    # m_values = [0]
    m_converged = []

    # if config.misc['parallel']:
    #     with multiprocessing.Pool() as p:
    #         images = p.map(compute_moment1, m_values)
    # else:
    #     for m_val in m_values:
    #         images.append(compute_moment1(m_val))
    for m_it in m_init:
        error = 1
        while abs(error) > 1e-4:
            m_new = compute_moment1(m_it)
            error = m_new - m_it
            m_it = m_new
        m_converged.append(m_it)
    return m_converged

    # fig, ax = plt.subplots(1, 1)
    # ax.plot(m_values, m_values, label='m')
    # ax.plot(m_values, images, label='R(m)')
    # ax.legend()
    # plt.show()


def convergence():
    vprint("[convergence]")

    factor_eval = quad_num.discretize(factor)

    def get_ground_state(eig_vec):
        ground_series = eig_vec * np.sign(eig_vec.coeffs[0])
        ground_state_eval = quad_num.eval(ground_series)
        ground_state_eval = ground_state_eval * factor_eval
        norm = quad_num.integrate(ground_state_eval, l2=True)
        ground_state_eval = ground_state_eval / norm
        return ground_state_eval

    vprint("-- Calculating exact solution")
    if params['Vp'].value.diff(x, x, x) == 0 and params['θ'].value == 0:
        solution = eq.solve_gaussian(forward, f, [x, y])
        solution_eval = quad_num.discretize(solution)
    else:
        eig_vals, eig_vecs = r_mat.eigs(k=1, which='LR')
        solution_eval = get_ground_state(eig_vecs[0])
    norm_sol = quad_num.integrate(solution_eval, l2=True)
    assert abs(norm_sol - 1) < 1e-6

    v0 = None
    degrees = []
    errors = []
    for d in range(5, degree):
        print("-- Solving for degree = " + str(d))
        npolys = core.iterator_size(2, d)
        sub_varf = r_mat.subdegree(d)
        # sub_mat = sub_mat
        if v0 is not None:
            actual_v0 = np.zeros(npolys)
            for i in range(len(v0)):
                actual_v0[i] = v0[i]
            v0 = actual_v0
        eig_vals, eig_vecs = sub_varf.eigs(v0=v0, k=1, which='LR')
        v0 = eig_vecs[0].coeffs.copy(order='C')
        ground_state_eval = get_ground_state(eig_vecs[0])
        error = la.norm(ground_state_eval - solution_eval, 2)
        degrees.append(d)
        errors.append(error)

    fig, ax = plt.subplots(1, 1)
    ground_series = eig_vecs[0] * np.sign(eig_vecs[0].coeffs[0])
    cont = quad_visu.plot(ground_series, factor, ax=ax)
    plt.colorbar(cont, ax=ax)
    plt.show()
    fig, ax = plt.subplots(1, 1)
    ax.semilogy(degrees, errors, 'k.')
    plt.show()
    # splot.contour(solution, (x, -1, 1), (y, -1, 1))


convergence()

factor_sym = factor
delta_arc_length = 0.1
beta = 10
# beta_range = np.arange(10, 1, delta_beta)
betas, conv, conv1, conv2 = [], [], [], []
m_init = [-2, 2]

while beta > 1:
    betas.append(beta)
    print("\n--> β = " + str(beta))
    backward_b = backward.subs(params['β'].symbol, beta).simplify()

    # vprint("Splitting operator in m-linear and m-independent parts")
    m_operator = backward_b.diff(params['m'].symbol)
    r_operator = (backward_b - params['m'].symbol*m_operator).cancel()

    # vprint("Discretizing operators")
    m_mat = quad_num.discretize_op(m_operator, f, degree, 2)
    r_mat = quad_num.discretize_op(r_operator, f, degree, 2)
    eig_vals, eig_vecs = r_mat.eigs(k=1, which='LR')

    bif = hermipy.cache()(bifurcation)
    conv.append(bifurcation(factor_x.subs(params['β'].symbol, beta), m_init))
    if len(conv) > 1:
        for i in range(len(m_init)):
            length = np.sqrt((betas[-1] - betas[-2])**2 + (conv[-1][i] - conv[-2][i])**2)
            delta_beta = (betas[-1] - betas[-2]) * delta_arc_length / length
            slope = (conv[-1][i] - conv[-2][i]) / (betas[-1] - betas[-2])
            beta = betas[-1] + delta_beta
            m_init[i] = conv[-1][i] + slope * delta_beta
    else:
        m_init = conv[-1]
        beta = beta - delta_arc_length
    # print(conv)
    # print(m_init)
fig, ax = plt.subplots(1, 1)
conv1 = [c[0] for c in conv]
conv2 = [c[1] for c in conv]
ax.plot(betas, conv1, 'k.')
ax.plot(betas, conv2, 'k.')
plt.show()

import ipdb;
ipdb.set_trace()
