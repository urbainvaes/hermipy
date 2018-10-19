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

# IMPORT MODULES {{{
import pdb
import argparse
import math
import os
import sympy as sym
import importlib
import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as las

import hermipy
import hermipy.stats
import hermipy.equations as eq
import hermipy.core as core

import matplotlib
from sympy.printing import pprint

sym.init_printing()

# }}}
# PARSE OPTIONS {{{
parser = argparse.ArgumentParser()

parser.add_argument('-tc', '--convergence', action='store_true',
                    help='Run convergence test')
parser.add_argument('-tce', '--convergence_eig', action='store_true',
                    help='Run convergence test with eig')
parser.add_argument('-tb', '--bifurcation', action='store_true',
                    help='Run bifurcation test')
parser.add_argument('-tp', '--plots', action='store_true')
parser.add_argument('-tw', '--white', action='store_true',
                    help='Generate bifurcation diagram for white noise')
parser.add_argument('-te0', '--epsilon0', action='store_true',
                    help='Study ε --> 0')

parser.add_argument('-v', '--verbose', action='store_true',
                    help='Enable verbose output')
parser.add_argument('--cache', action='store_true',
                    help='Enable use of cache')
parser.add_argument('-i', '--interactive', action='store_true',
                    help='Enable interactive plots')

parser.add_argument('-c', '--config', type=str, help='Configuration file')
parser.add_argument('-b', '--beta', type=str, help='Value of β')
parser.add_argument('-t', '--theta', type=str, help='Value of θ')
parser.add_argument('-e', '--epsilon', type=str, help='Value of ε')
parser.add_argument('-g', '--gamma', type=str, help='Value of γ')
parser.add_argument('-m', '--mass', type=str, help='Value of m')
parser.add_argument('-m0', '--m0', type=str, help='Initial value of m')
parser.add_argument('-b0', '--b0', type=str, help='Initial value of β')
parser.add_argument('-d', '--directory', type=str, help='Directory of output')

args = parser.parse_args()

if args.config:
    config = importlib.import_module(args.config)
else:
    config = importlib.import_module("examples.config")

# Make sure contains arrays
try:
    config.eq
except NameError:
    config.eq = {}

try:
    config.num
except NameError:
    config.num = {}

try:
    config.misc
except NameError:
    config.misc = {}

if args.directory:
    config.misc['directory'] = args.directory + '/'
if 'directory' not in config.misc:
    config.misc['directory'] = './'
if not os.path.exists(config.misc['directory']):
    os.makedirs(config.misc['directory'])

if args.beta:
    config.eq['β'] = sym.Rational(args.beta)

if args.theta:
    config.eq['θ'] = sym.Rational(args.theta)

if args.epsilon:
    config.eq['ε'] = sym.Rational(args.epsilon)

if args.gamma:
    config.eq['γ'] = sym.Rational(args.gamma)

if args.mass:
    config.eq['m'] = sym.Rational(args.mass)

if args.m0:
    config.num['m0'] = float(args.m0)

if args.verbose:
    config.misc['verbose'] = args.verbose
elif 'verbose' not in config.misc:
    config.misc['verbose'] = False

if not args.interactive:
    matplotlib.use('Agg')

if True:  # To placate linter
    import matplotlib.pyplot as plt

if args.cache:
    config.misc['cache'] = args.cache
elif 'cache' not in config.misc:
    config.misc['cache'] = False


def vprint(*args, **kwargs):
    if config.misc['verbose']:
        print(*args, **kwargs)


# Set library option
hermipy.settings.update(config.misc)
vprint(hermipy.settings)

matplotlib.rc('font', size=14)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)

# }}}
# DATA AND PARAMETERS FOR NUMERICAL SIMULATION {{{
vprint("Symbolic calculations")

# Variables and function
x, y = eq.McKean_Vlasov.x, eq.McKean_Vlasov.y
f = eq.McKean_Vlasov.f

# Parameters
params = {}


class Parameter():

    def __init__(self, symbol, value, type="main"):
        self.symbol = symbol
        self.value = value
        self.type = type

    def __repr__(self):
        return "Symbol: " + str(self.symbol) + " - " \
               + "Value: " + str(self.value) + "."

    def eval(self):
        value = self.value
        sym_val = {params[k].symbol: params[k].value for k in params}
        m, β = params['m'].symbol, params['β'].symbol
        ε, γ = params['ε'].symbol, params['γ'].symbol
        while value.free_symbols - {x, y, m, β, ε, γ} != set():
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
    if 'Vy' in config.eq:
        eq_params['Vy'] = Parameter(symbol, config.eq['Vy'], type='equation')
    return eq_params


def vy():
    r, symb = sym.Rational, sym.symbols
    if 'Vy' in config.eq:
        r, symb, β = sym.Rational, sym.symbols, params['β']
        μy = Parameter(symb('μy', real=True), config.num['μy'])
        σy = Parameter(symb('σy', real=True, positive=True), config.num['σy'])
        Vqy = Parameter(sym.Function('Vqy')(y), r(1, 2) *
                        (y - μy.symbol)*(y - μy.symbol)/(σy.symbol))
        return {'μy': μy, 'σy': σy, 'Vqy': Vqy}
    else:
        Vy = Parameter(sym.Function('Vy')(y), y*y/2, type='equation')
        return {'Vqy': Vy, 'Vy': Vy}


def vq():
    r, symb, β = sym.Rational, sym.symbols, params['β']
    μx = Parameter(symb('μx', real=True), config.num['μx'])
    σx = Parameter(symb('σx', real=True, positive=True), config.num['σx'])
    Vq = Parameter(sym.Function('Vqx')(x), r(1, 2) *
                   (x - μx.symbol)*(x - μx.symbol)/(σx.symbol))
    return {'μx': μx, 'σx': σx, 'Vqx': Vq}


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
    forward = eq.McKean_Vlasov.equation(eq_params)
    fluxes = eq.McKean_Vlasov.fluxes(eq_params)
    return forward, fluxes


def factors(symbolic, λ):
    β = params['β'].value if symbolic == 0 else params['β'].symbol
    if symbolic == 2:
        Vp, Vq = params['Vp'].symbol, params['Vqx'].symbol
        Vy, Vqy = params['Vy'].symbol, params['Vqy'].symbol
    elif symbolic == 1:
        Vp, Vq = params['Vp'].value, params['Vqx'].value
        Vy, Vqy = params['Vy'].value, params['Vqy'].value
    elif symbolic == 0:
        Vp, Vq = params['Vp'].eval(), params['Vqx'].eval()
        Vy, Vqy = params['Vy'].eval(), params['Vqy'].eval()
    if args.bifurcation:
        factor_x = sym.exp(-Vq/2)
        factor_y = sym.exp(-Vqy/2)
    else:
        factor_x = sym.exp(-(λ*Vq + β*(1-λ)*Vp))
        factor_y = sym.exp(-(λ*Vqy + (1-λ)*Vy))
    factor_x = sym.exp(-Vq/2)
    factor_y = sym.exp(-(λ*Vqy + (1-λ)*Vy))
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
forward, fluxes = forward_op(config.misc['symbolic'])

# Mapped backward operator
factor_x, factor_y, factor = factors(config.misc['symbolic'], config.num['λ'])
backward = eq.map_operator(forward, f, factor)

forward, backward = evaluate([forward, backward])
fluxes = evaluate(fluxes)
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
    μy = params['μy'].value
    σx = params['σx'].value
    σy = params['σy'].value

    quad_num = hermipy.Quad.gauss_hermite(
            config.num['n_points_num'], dim=2,
            mean=[μx, μy], cov=[[σx, 0], [0, σy]],
            factor=factor)

    # For visualization
    nv = 200
    # band_width = np.sqrt(2) * np.sqrt(2*degree + 1)
    # bounds_x = band_width*np.sqrt(float(σx))
    # bounds_y = band_width*np.sqrt(float(σy))
    bounds_x = 4
    bounds_y = 4

    quad_visu = hermipy.Quad.newton_cotes([nv, nv], [bounds_x, bounds_y],
                                          factor=factor)
    return quad_num, quad_visu


quad_num, quad_visu = compute_quads()

# }}}
# SPECTRAL METHOD FOR STATIONARY EQUATION {{{
vprint("Splitting operator in m-linear and m-independent parts")
index_set = config.num['index_set']
m_operator = forward.diff(params['m'].symbol)
r_operator = (forward - params['m'].symbol*m_operator).cancel()
if params['m'].value.free_symbols == set() \
        and params['β'].value.free_symbols == set() \
        and params['ε'].value.free_symbols == set() \
        and params['γ'].value.free_symbols == set():
    fd = [quad_num.discretize_op(flux, degree, index_set=index_set) for flux in fluxes]

vprint("Solving eigenvalue problem")

# Make plots {{{


def plot():

    fd = [quad_num.discretize_op(flux, degree, index_set=index_set) for flux in fluxes]
    r_mat = quad_num.discretize_op(r_operator, degree, index_set=index_set)
    eig_vals, eig_vecs = r_mat.eigs(k=1, which='LR')
    vprint("[plot]")

    fig, ax = plt.subplots(1, 1)
    ground_series = eig_vecs[0] * np.sign(eig_vecs[0].coeffs[0])
    ground_series = ground_series / quad_num.integrate(ground_series, flat=True)
    cont = quad_visu.plot(ground_series, ax=ax, vmin=0,
                          extend='min', bounds=False)
    for c in cont.collections:
        c.set_edgecolor("face")
    plt.colorbar(cont, ax=ax)
    f0, f1 = fd[0](ground_series), fd[1](ground_series)
    streamlines = quad_visu.streamlines(f0, f1, ax=ax)
    for c in streamlines.collections:
        plt.setp(c, linewidth=.5)
    output = config.misc['directory'] + 'solution.eps'
    plt.savefig(output, bbox_inches='tight')
    plt.show()
    fig, ax = plt.subplots(1, 1)
    scatter = ground_series.subdegree(10).plot(ax=ax)
    plt.colorbar(scatter, ax=ax)
    output = config.misc['directory'] + 'coefficients.eps'
    plt.savefig(output, bbox_inches='tight')
    plt.show()

    def plot_eigenfunctions():
        fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
        axes = (ax11, ax12, ax21, ax22)
        for e_val, series, ax in zip(eig_vals, eig_vecs, axes):
            series = series * np.sign(series.coeffs[0])
            kwargs = {} if abs(e_val) > 1e-5 else {'vmin': 0, 'extend': 'min'}
            cont = quad_visu.plot(series, ax=ax, bounds=False, **kwargs)
            ax.set_title("Eigenvalue: " + str(e_val))
            plt.colorbar(cont, ax=ax)
        plt.show()

    def plot_hermite_functions():
        quad_num.plot_hf([0, 0])
        quad_num.plot_hf([degree, 0])
        quad_num.plot_hf([0, degree])
        quad_num.plot_hf([degree // 2, degree // 2])

    def plot_ground_state():
        ground_state = eig_vecs[0] * np.sign(eig_vecs[0].coeffs[0])

        # Plot the projections
        degrees = np.arange(degree + 1)
        series_x = ground_state.project(0)
        series_y = ground_state.project(1)
        quad_x = quad_visu.project(0)
        quad_y = quad_visu.project(1)
        fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
        quad_x.plot(series_x, ax=ax11)
        ax11.set_title("Marginal distribution in x direction")
        ax12.bar(degrees, series_x.coeffs)
        ax12.set_title("Hermite coefficients of x marginal")
        quad_y.plot(series_y, ax=ax21)
        ax21.set_title("Marginal distribution in y direction")
        ax22.bar(degrees, series_y.coeffs)
        ax22.set_title("Hermite coefficients of y marginal")
        plt.show()

    def plot_comparison_with_asym():
        fig, ax = plt.subplots(1, 1)
        as_sol = sym.exp(- params['β'].value * params['Vp'].eval())
        as_series = quad_num.project(0).transform(as_sol, degree)
        quad_visu.project(0).plot(as_series, ax=ax)
        plt.show()

    # def plot_discretization_error():
    #     fig, ax = plt.subplots(1, 2)
    #     quad_visu_x = quad_visu.project(0)
    #     as_sol = sym.exp(- params['β'].value * params['Vp'].eval())
    #     as_series = quad_num.project(0).transform(as_sol/factor_x, degree)
    #     norm_as = la.norm(as_series.coeffs, 2)
    #     as_series.coeffs = as_series.coeffs / norm_as
    #     quad_visu_x.plot(as_series, degree, factor_x, ax[0])
    #     sol_x = quad_visu_x.discretize(as_sol)
    #     ax[0].plot(quad_visu_x.discretize('x'), sol_x/norm_as)
    #     plt.show()

    plot_eigenfunctions()
    plot_hermite_functions()
    plot_ground_state()
    plot_comparison_with_asym()


if args.plots:
    plot()
# }}}
# Test convergence {{{

def convergence_eig():
    vprint("[convergence]")
    vprint("Discretizing operators")
    r_mat = quad_num.discretize_op(r_operator, degree, index_set=index_set)
    # r_mat.plot()

    def asymptotic():
        Vp = params['θ'].value*(x - params['m'].value)**2/2
        Vp = Vp + params['Vp'].eval()
        β, ε = params['β'].value, params['ε'].value
        rho_eta = 1/sym.sqrt(2*np.pi) * sym.exp(-y*y/2)
        quadx = quad_num.project(0)
        Z = quadx.integrate(sym.exp(-β*Vp), flat=True)
        u0 = sym.exp(-β*Vp)/Z
        u1 = y*sym.sqrt(β)*sym.exp(-β*Vp)*Vp.diff(x)/Z
        C2 = quadx.integrate(
                (Vp.diff(x)**2/2 - 1/β*Vp.diff(x, x))
                * sym.exp(-β*Vp)/Z, flat=True)
        u2 = (2*C2*β + y**2*β*Vp.diff(x)**2 - y**2*Vp.diff(x, x)
              - 2*β*Vp.diff(x)**2 + 3*Vp.diff(x, x))*sym.exp(-β*Vp)/(2*Z)
        u3 = u0.diff(x, x, x)/(6*float(β)**(3/2))*(-y**3 - 3*y) \
            - y/np.sqrt(float(β))*(Vp.diff(x)*u0.diff(x)).diff(x) \
            - y/np.sqrt(float(β))*u2.diff(x)
        ua = (u0 + ε*u1 + ε**2*u2)*rho_eta
        # ta = quad_num.transform(ua/factor, degree=degree)
        # quad_visu.plot(ua)
        assert abs(quad_num.integrate(ua, flat=True) - 1) < 1e-6
        return quad_num.discretize(ua)

    asym_num = asymptotic()

    def get_ground_state(eig_vec):
        ground_series = eig_vec * np.sign(eig_vec.coeffs[0])
        # quad_visu.plot(ground_series, factor=factor)
        ground_state_eval = quad_num.eval(ground_series)
        norm = quad_num.norm(ground_state_eval, n=1, flat=True)
        ground_state_eval = ground_state_eval / norm
        return ground_state_eval

    vprint("-- Calculating exact solution")
    if params['Vp'].value.diff(x, x, x) == 0 and params['θ'].value == 0:
        solution = eq.solve_gaussian(forward, f, [x, y])
        solution_eval = quad_num.discretize(solution)
        norm_sol = quad_num.integrate(solution_eval, flat=True)
        assert abs(norm_sol - 1) < 1e-6
    else:
        eig_vals, eig_vecs = r_mat.eigs(k=4, which='LR')
        argmin = np.argmin(np.abs(eig_vals))
        solution_eval = get_ground_state(eig_vecs[argmin])
        norm_sol = quad_num.integrate(solution_eval, flat=True)
        solution_eval = solution_eval / norm_sol

    quad_num.plot(solution_eval)

    v0, degrees, errors, errors_a, mins, eig_ground = None, [], [], [], [], []
    for d in range(20, degree):
        print("-- Solving for degree = " + str(d))
        npolys = core.iterator_size(2, d, index_set=index_set)
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
        # if d > 60:
        #     quad_num.plot(eig_vecs[0], title=str(eig_vals[0]))
        # error_consistency =
        # error_discretization =
        error = quad_num.norm(ground_state_eval - solution_eval, n=1, flat=True)
        error_a = quad_num.norm(ground_state_eval - asym_num, n=1, flat=True)
        degrees.append(d)
        mins.append(np.min(ground_state_eval))
        eig_ground.append(eig_vals[0])
        errors.append(error)
        errors_a.append(error_a)
        print(error, eig_vals[0], np.min(ground_state_eval))

    fig, ax = plt.subplots(1, 1)
    ground_series = eig_vecs[0] * np.sign(eig_vecs[0].coeffs[0])
    cont = quad_visu.plot(ground_series, ax=ax,
                          vmin=0, extend='min', bounds=False)
    plt.colorbar(cont, ax=ax)
    plt.show()

    def plot_log(x, y, file, lin=True):
        x, y = np.asarray(x), np.asarray(y)
        fig, ax = plt.subplots(1, 1)
        ax.semilogy(x, y, 'k.', label='$|\\lambda_0|$')
        ax.set_yscale('log', basey=2)
        cut_off = 70
        x_poly = x[0:cut_off + 1] if len(x) > cut_off else x
        y_poly = y[0:cut_off + 1] if len(y) > cut_off else y
        coeffs = np.polyfit(x_poly, np.log10(y_poly), 1)
        if lin:
            ax.semilogy(x, 10**coeffs[1] * 10**(coeffs[0]*x), 'r-',
                        label='$y = {:.2f} \\times 10^{{ {:.2f} \\, x}}$'.
                              format(10**coeffs[1], coeffs[0]))
        plt.legend(loc='upper right')
        ax.set_xlabel('d')
        plt.savefig(file, bbox_inches='tight')
        plt.show()

    # np.save('errors-epsilon-'+str(float(params['ε'].value)), errors_a)
    degrees = degrees[0:-10]
    mins = mins[0:-10]
    errors = errors[0:-10]
    eig_ground = eig_ground[0:-10]

    dir = config.misc['directory']
    plot_log(degrees, errors, dir + 'l1conv.eps')
    plot_log(degrees, [max(0, -m) for m in mins], dir + 'min.eps')
    plot_log(degrees, np.abs(eig_ground), dir + 'eig_val.eps', lin=False)
    # plot_log(degrees, errors_a, dir + 'l1asym.eps', lin=False)

    # qy = quad_num.project(1)
    # wy = qy.factor * qy.factor / qy.position.weight()
    # Iy = qy.transform(wy, degree=degree, index_set=index_set)
    # quad_num.plot(Iy*eig_vecs[0])

    fig, ax = plt.subplots()
    ax.semilogy(degrees, errors, '.', label="$L^1$ error")
    ax.semilogy(degrees, [max(0, -m) for m in mins], '.', label="- minimum")
    # ax.semilogy(degrees, np.abs(eig_ground), '.', label="$|\\lambda_0|$")
    # condition = (np.asarray(degrees) + 1) % 2
    # xplot, yplot = np.extract(condition, degrees), np.extract(condition, errors)
    # ax.semilogy(xplot, yplot, 'b.', label="$L^1$ error")
    # coeffs = np.polyfit(xplot, np.log10(yplot), 1)
    # ax.semilogy(xplot, 10**coeffs[1] * 10**(coeffs[0]*xplot), 'b-')
    # yplot = np.extract(condition, np.abs(eig_ground))
    # ax.semilogy(xplot, yplot, 'r.', label="$|\\lambda_0(d)|$")
    # coeffs = np.polyfit(xplot, np.log10(yplot), 1)
    # ax.semilogy(xplot, 10**coeffs[1] * 10**(coeffs[0]*xplot), 'r-')
    plt.legend(loc='upper right')
    ax.set_xlabel('d')
    plt.savefig("errors.eps", bbox_inches='tight')
    plt.show()


if args.convergence_eig:
    convergence_eig()
# }}}
# {{{ Time dependent solution

# Correction to the drift
if 'drift_correction' in config.num:
    drift = config.num['drift_correction']
    scaling = sym.sqrt(2/params['β'].value)/params['ε'].value
    forward = forward - scaling * drift * f.diff(x)


def white_noise_bifurc():
    βmin, βmax, sstep = 1, 15, .01
    qx = quad_num.project(0)

    def branch(m0):
        β, m, βs, ms = βmax, m0, [], []
        Vx, θ = params['Vp'].eval(), params['θ'].value
        while β > βmin:
            threshold, Δ, Δold, mold = 1e-5, None, None, None
            while Δ is None or abs(Δ) > threshold:
                Vt = θ*(x - m)**2/2 + Vx
                u = sym.exp(-β*Vt)
                integral = qx.integrate(u, flat=True)
                Δ = qx.integrate(u/integral*hermipy.x, flat=True) - m
                if Δold is None:
                    newm = m + Δ
                else:
                    # secant's method
                    newm = (mold*Δ - m*Δold)/(Δ - Δold)
                m, mold, Δold = newm, m, Δ
                print(β, m)
            βs.append(β)
            ms.append(m)
            gmm, dsdβ = 20, 1
            if len(ms) > 1:
                Δm, Δβ = ms[-1] - ms[-2], βs[-1] - βs[-2]
                dsdβ = np.sqrt(gmm*Δm*Δm + Δβ*Δβ) / abs(Δβ)
            β = β - sstep/dsdβ
        return ms, βs

    m_up, b_up = branch(1)
    m_down, b_down = branch(-1)
    dir = config.misc['directory']
    np.save(dir + "white-noise-ms-up", np.asarray(m_up))
    np.save(dir + "white-noise-betas-up.npy", np.asarray(b_up))
    np.save(dir + "white-noise-ms-down.npy", np.asarray(m_down))
    np.save(dir + "white-noise-betas-down.npy", np.asarray(b_down))


if args.white:
    white_noise_bifurc()


def time_dependent():

    use_translate = False
    t_sym, Vp = sym.symbols('t'), params['Vp'].value
    translate_Vp = (Vp.subs(x, x + t_sym) - Vp).expand()
    translate_op = (translate_Vp.diff(x)*f).diff(x)

    # Output directory
    dir = config.misc['directory']

    if args.interactive:
        plt.ion()
        fig, ax = plt.subplots(2, 2)

    m_operator = forward.diff(params['m'].symbol)
    r_operator = (forward - params['m'].symbol*m_operator).cancel()
    m_mat = quad_num.discretize_op(m_operator, degree, index_set=index_set)
    sstep, βmin, βmax = .05, 0.4, sym.Rational(args.b0) if args.b0 else 10

    # Calculate projections
    qx, qy = quad_num.project(0), quad_num.project(1)
    wx = qx.factor * qx.factor / qx.position.weight()
    wy = qy.factor * qy.factor / qy.position.weight()

    # Integral and moment operators
    Ix = qx.transform(wx, degree=degree, index_set=index_set)
    Iy = qy.transform(wy, degree=degree, index_set=index_set)
    mx1 = qx.transform(wx * x, degree=degree, index_set=index_set)
    my1 = qy.transform(wy * y, degree=degree, index_set=index_set)

    assert 'm0' in config.num
    β, m, betas, ms = βmax, config.num['m0'], [], []
    for i in range(20):
        print("m = {}".format(m))
        Vp = params['θ'].value*(x - m)**2/2
        Vp = Vp + params['Vp'].eval()
        u = sym.exp(-β*Vp)
        integral = qx.integrate(u, flat=True)
        m = qx.integrate(u / integral * hermipy.x, flat=True)
    m = round(m, 2)

    if use_translate:
        translation = m
        u = u.subs(x, x + translation)
        m = 0
    else:
        translation = 0

    u = u * sym.exp(-params['Vy'].eval())
    integral = quad_num.integrate(u, flat=True)
    t = quad_num.transform(u/integral, degree=degree, index_set=index_set)

    x_eval = quad_num.discretize('x')
    eye = quad_num.varf('1', degree=degree, index_set=index_set)
    dt, Ns, scheme = 2**-9, int(2e3), "backward"

    while β > βmin:

        kwargs = {'degree': degree, 'index_set': index_set}
        r_operator_this = r_operator.subs(params['β'].symbol, β)
        r_mat = quad_num.discretize_op(r_operator_this, **kwargs)

        if use_translate:
            t_operator_this = translate_op.subs(t_sym, translation)
            r_mat = r_mat + quad_num.discretize_op(t_operator_this, **kwargs)

        Vx = params['Vp'].eval().subs(x, x + translation)
        θ = params['θ'].value

        difference = np.inf
        for i in range(Ns):

            if i % 10 == 0 and args.interactive:

                ax[0][0].clear()
                quad_visu.plot(t, ax=ax[0][0], bounds=True,
                               vmin=0, extend='min')

                ax[0][1].clear()
                density = sym.exp(-β*(Vx + θ*(x - m)**2/2))
                density = density / qx.integrate(density, flat=True)

                t.plot(ax=ax[0][1])

                # Projections
                proj_x, proj_y = Iy*t, Ix*t
                mx, my = mx1*proj_x, my1*proj_y
                ax[1][0].clear()
                quad_visu.project(0).plot(proj_x, ax=ax[1][0])
                quad_visu.project(0).plot(density, ax=ax[1][0])
                ax[1][0].set_title("First moment: " + str(mx.coeffs[0]))

                ax[1][1].clear()
                quad_visu.project(1).plot(proj_y, ax=ax[1][1])
                ax[1][1].set_title("First moment: " + str(my.coeffs[0]))

                plt.draw()
                plt.pause(.01)

            print("β: " + str(β) + ", i: " + str(i), "dt: " + str(dt) +
                  ", m: " + str(m + translation), ", Δ: " + str(difference))

            operator = r_mat + m*m_mat

            # Backward Euler
            if scheme == 'backward':
                total_op = eye - dt*operator
                new_t = total_op.solve(t)

            # Crank-Nicholson
            if scheme == 'crank':
                crank_left = eye - dt*operator*(1/2)
                crank_right = eye + dt*operator*(1/2)
                new_t = crank_left.solve(crank_right(t))

            # Normalization
            r_eval = quad_num.eval(new_t)
            integral = quad_num.integrate(r_eval, flat=True)
            new_t, r_eval = new_t * (1/integral), r_eval * (1/integral)
            new_m = quad_num.integrate(r_eval*x_eval, flat=True)

            difference = quad_num.norm(t - new_t, n=1)/dt
            difference = difference + abs(new_m - m)/dt

            # Time adaptation
            threshold, dt_max = .01, 256

            if difference*dt > 2*threshold:
                dt = dt / 2.
                continue
            elif difference*dt < threshold and dt < dt_max:
                dt = dt * 2.
            t, m = new_t, new_m

            if difference < 5e-8:
                betas.append(β)
                ms.append(m + translation)

                gmm, dsdβ = 20, 1
                if len(ms) > 1:
                    Δm, Δβ = ms[-1] - ms[-2], betas[-1] - betas[-2]
                    dsdβ = np.sqrt(gmm*Δm*Δm + Δβ*Δβ) / abs(Δβ)
                newβ = β - sstep/dsdβ

                if math.floor(β) != math.floor(newβ):
                    β = math.floor(β)

                    fig, ax = plt.subplots(1, 1)
                    quad_visu.plot(t, bounds=False, ax=ax, title="$\\rho(x, \\eta)$")
                    plt.savefig(dir + 'solution-beta=' + str(β) + '.eps',
                                bbox_inches='tight')
                    plt.close()

                    fig, ax = plt.subplots(1, 1)
                    density = sym.exp(-β*(Vx + θ*(x - m)**2/2))
                    density = density / qx.integrate(density, flat=True)
                    quad_visu.project(0).plot(density, ax=ax, label="White noise")
                    ε = params['ε'].value
                    quad_visu.project(0).plot(Iy*t, bounds=False, ax=ax, label="$\\varepsilon = " + str(ε) + "$")
                    plt.legend()
                    ax.set_title("$\\int \\rho(x, \\eta) \\, \\mathrm d \\eta$")
                    plt.savefig(dir + 'solution-proj-beta=' + str(β) + '.eps',
                                bbox_inches='tight')
                    plt.close()

                β = newβ

                if use_translate:
                    translation += m
                    m = 0

                break

    ε = config.eq['ε']
    np.save(dir + "epsilon=" + str(ε).replace('/', 'o') + "-m0=" + str(args.m0) + "-betas", np.asarray(betas))
    np.save(dir + "epsilon=" + str(ε).replace('/', 'o') + "-m0=" + str(args.m0) + "-ms", np.asarray(ms))


if args.bifurcation:
    time_dependent()


def convergence():

    # Output directory
    dir = config.misc['directory']

    if args.interactive:
        plt.ion()
        fig, ax = plt.subplots(2, 2)

    r_operator = (forward - params['m'].symbol*m_operator).cancel()
    r_mat = quad_num.discretize_op(r_operator, degree, index_set=index_set)

    # Calculate projections
    qx, qy = quad_num.project(0), quad_num.project(1)
    wx = qx.factor * qx.factor / qx.position.weight()
    wy = qy.factor * qy.factor / qy.position.weight()

    # Integral and moment operators
    Ix = qx.transform(wx, degree=degree, index_set=index_set)
    Iy = qy.transform(wy, degree=degree, index_set=index_set)

    m, Vx = params['m'].value, params['Vp'].eval()
    β, θ = params['β'].value, params['θ'].value
    Vp = Vx + θ*(x - m)**2/2
    u = sym.exp(-β*Vp)*sym.exp(-params['Vy'].eval())
    integral = quad_num.integrate(u, flat=True)
    t = quad_num.transform(u/integral, degree=degree, index_set=index_set)

    eye = quad_num.varf('1', degree=degree, index_set=index_set)
    dt, Ns, scheme = 2**-9, int(1e4), "backward"

    density = sym.exp(-β*(Vx + θ*(x - m)**2/2))
    density = density / qx.integrate(density, flat=True)

    dmin, d, degrees, errors = 20, degree, [], []
    while d >= dmin:
        r_mat, eye = r_mat.subdegree(d), eye.subdegree(d)
        t, Ix, Iy = t.subdegree(d), Ix.subdegree(d), Iy.subdegree(d)

        Δ = np.inf
        for i in range(Ns):

            print("d: " + str(d) + ", i: " + str(i) +
                  ", dt: " + str(dt) + ", Δ: " + str(Δ))

            # Backward Euler
            if scheme == 'backward':
                total_op = eye - dt*r_mat
                new_t = total_op.solve(t)

            # Crank-Nicholson
            if scheme == 'crank':
                crank_left = eye - dt*r_mat*(1/2)
                crank_right = eye + dt*r_mat*(1/2)
                new_t = crank_left.solve(crank_right(t))

            # Normalization
            r_eval = quad_num.eval(new_t)
            integral = quad_num.integrate(r_eval, flat=True)
            new_t = new_t * (1/integral)

            Δ = quad_num.norm(t - new_t, n=1)/dt

            # Time adaptation
            threshold, dt_max = .01, 64
            if Δ*dt > 2*threshold:
                dt = dt / 2.
                continue
            elif Δ*dt < threshold and dt < dt_max:
                dt = dt * 2.
            t = new_t

            if Δ < 1e-10:

                if d == config.num['degree']:
                    t_exact = t
                    e_exact = quad_num.eval(t_exact)

                else:
                    degrees.append(d)
                    error = quad_num.eval(t) - e_exact
                    errors.append(quad_num.norm(error, n=1, flat=True))
                    print(errors[-1])

                if args.interactive:

                    ax[0][0].clear()
                    quad_visu.plot(t, ax=ax[0][0], bounds=True,
                                   vmin=0, extend='min')

                    # Hermite coefficients
                    ax[0][1].clear()
                    t.plot(ax=ax[0][1])

                    # Projections
                    proj_x, proj_y = Iy*t, Ix*t
                    ax[1][0].clear()
                    quad_visu.project(0).plot(proj_x, ax=ax[1][0])
                    quad_visu.project(0).plot(density, ax=ax[1][0])

                    ax[1][1].clear()
                    quad_visu.project(1).plot(proj_y, ax=ax[1][1])

                    plt.draw()
                    plt.pause(.01)

                # fig, ax = plt.subplots(1)
                # cont = quad_visu.plot(t, bounds=False, ax=ax,
                #                       title="$\\rho(x, \\eta)$")
                # plt.colorbar(cont, ax=ax)
                # plt.savefig('solution.eps', bbox_inches='tight')
                # plt.close()

                # fig, ax = plt.subplots(1)
                # title = "Hermite coefficients"
                # hermite_coeffs = t.plot(ax=ax, title=title)
                # plt.colorbar(hermite_coeffs, ax=ax)
                # plt.savefig('hermite_coefficients.eps', bbox_inches='tight')
                # plt.close()

                d = d - 1
                break
    import ipdb; ipdb.set_trace()


if args.convergence:
    convergence()


def convergence_epsilon():

    # Output directory
    dir = config.misc['directory']

    if args.interactive:
        plt.ion()
        fig, ax = plt.subplots(2, 2)

    # Calculate projections
    qx, qy = quad_num.project(0), quad_num.project(1)
    wx = qx.factor * qx.factor / qx.position.weight()
    wy = qy.factor * qy.factor / qy.position.weight()

    # Integral and moment operators
    Ix = qx.transform(wx, degree=degree, index_set=index_set)
    Iy = qy.transform(wy, degree=degree, index_set=index_set)
    mx1 = qx.transform(wx * x, degree=degree, index_set=index_set)
    my1 = qy.transform(wy * y, degree=degree, index_set=index_set)

    β = params['β'].eval()
    Vx, Vy = params['Vp'].eval(), params['Vy'].eval()
    ux, uy = sym.exp(-β*Vx), sym.exp(-Vy)
    Zx, Zy = qx.integrate(ux, flat=True), qy.integrate(uy, flat=True)
    ux, uy, u = ux/Zx, uy/Zy, ux/Zx*uy/Zy
    integral = quad_num.integrate(u, flat=True)
    t = quad_num.transform(u/integral, degree=degree, index_set=index_set)

    eye = quad_num.varf('1', degree=degree, index_set=index_set)
    dt, Ns = 2**-9, int(1e4)

    # Calculate white noise solution
    #
    # Comparison with this is not great, because there is an error on the
    # coefficient that multiplies the noise.
    Vp, β = params['Vp'].value, params['β'].value
    γ, ε = params['γ'].symbol, params['ε'].symbol
    kwargs = {'degree': degree, 'index_set': index_set}
    forward0x = eq.Fokker_Planck_1d.equation({'Vp': Vp, 'β': β})
    var0x = qx.discretize_op(forward0x, **kwargs)
    _, [tx0] = var0x.eigs(k=1, which='LR')
    forward0 = forward.subs(γ, 1).subs(ε, 1)
    var0 = quad_num.discretize_op(forward0, **kwargs)
    _, [t20] = var0.eigs(k=1, which='LR')
    tx0 = tx0 * (1/qx.integrate(tx0, flat=True, tensorize=False))
    t20 = t20 * (1/quad_num.integrate(t20, flat=True, tensorize=False))

    εs = [2**(-i/4) for i in range(17)]
    e2, ex, t2, tx = [], [], [], []

    for iε, ε in enumerate(εs):

        kwargs = {'degree': degree, 'index_set': index_set}
        forward_ε = forward.subs(γ, 0).subs(params['ε'].symbol, ε)
        r_mat = quad_num.discretize_op(forward_ε, **kwargs)

        Δ, Δx = np.inf, np.inf
        for i in range(Ns):

            if i % 10 == 0 and args.interactive:

                ax[0][0].clear()
                quad_visu.plot(t, ax=ax[0][0], bounds=False,
                               vmin=0, extend='min')

                ax[0][1].clear()
                t.plot(ax=ax[0][1])

                # Projections
                proj_x, proj_y = Iy*t, Ix*t
                mx, my = mx1*proj_x, my1*proj_y
                ax[1][0].clear()
                quad_visu.project(0).plot(proj_x, ax=ax[1][0])
                quad_visu.project(0).plot(ux, ax=ax[1][0])
                ax[1][0].set_title("First moment: " + str(mx.coeffs[0]))

                ax[1][1].clear()
                quad_visu.project(1).plot(proj_y, ax=ax[1][1])
                quad_visu.project(1).plot(uy, ax=ax[1][1])
                ax[1][1].set_title("First moment: " + str(my.coeffs[0]))

                plt.draw()
                plt.pause(.01)

            print("ε: " + str(ε) + ", i: " + str(i), "dt: " + str(dt) +
                  ", Δ: " + str(Δ) + ", Δx: " + str(Δx))

            total_op = eye - dt*r_mat
            new_t = total_op.solve(t)

            # Normalization
            r_eval = quad_num.eval(new_t)
            integral = quad_num.integrate(r_eval, flat=True)
            new_t = new_t * (1/integral)

            Δ = quad_num.norm(t - new_t, n=1, flat=True)/dt
            Δx = qx.norm(Iy*(t - new_t), n=1, flat=True)/dt

            # Time adaptation
            threshold, dt_max = .01, 64

            if Δ*dt > 2*threshold:
                dt = dt / 2.
                continue
            elif Δ*dt < threshold and dt < dt_max:
                dt = dt * 2.

            t = new_t

            if Δ < 1e-16:
                t2.append(t)
                tx.append(Iy*t)
                error2, errorx = t20 - t2[-1], Iy*(t20 - t2[-1])
                e2.append(quad_num.norm(error2, n=1, flat=True))
                ex.append(qx.norm(errorx, n=1, flat=True))

                fig, ax = plt.subplots(1, 1)
                quad_visu.plot(t, ax=ax, vmin=0, extend='min')
                ax.set_title("$\\varepsilon = {}$".format(ε))
                plt.savefig(dir + 'solution-epsilon=' + str(ε) + '.eps',
                            bbox_inches='tight')
                plt.close()


                if iε > 2:
                    logε = np.log2(εs[0:iε+1])
                    coeffs = np.polyfit(logε, np.log2(e2), 1)
                    print(coeffs)
                    coeffs_x = np.polyfit(logε, np.log2(ex), 1)
                    print(coeffs_x)

                break

    fig, ax = plt.subplots()
    cmap = matplotlib.cm.get_cmap('viridis_r')
    for i, (ε, txε) in enumerate(zip(εs, tx)):
        if i % 2 is not 0:
            continue
        label = "$\\varepsilon = 2^{{ -{} }} $".format(i/4)
        kwargs = {'color': cmap(ε), 'label': label}
        quad_visu.project(0).plot(txε, ax=ax, **kwargs)
    ax.set_title("")
    plt.legend()
    plt.savefig("convergence-epsilon.eps", bbox_inches='tight')

    fig, ax = plt.subplots()
    xplot, yplot1 = logε = np.asarray(εs), np.asarray(e2)
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)
    ax.plot(xplot, yplot1, 'b.', label="$|\\rho - \\rho_0|_1$")
    yplot2 = np.asarray(ex)
    ax.plot(xplot, yplot2, 'r.', label="$|\\rho^x - \\rho^x_0|_1$")
    coeffs = np.polyfit(np.log2(xplot), np.log2(yplot1), 1)
    ax.plot(xplot, 2**coeffs[1] * xplot**coeffs[0], 'b-',
            label='$y = {:.2f} \\, \\times \\, \\varepsilon^{{ {:.2f} }}$'.
            format(2**coeffs[1], coeffs[0]))
    coeffs = np.polyfit(np.log2(xplot), np.log2(yplot2), 1)
    ax.plot(xplot, 2**coeffs[1] * xplot**coeffs[0], 'r-',
             label='$y = {:.2f} \\, \\times \\, \\varepsilon^{{ {:.2f} }}$'.
                              format(2**coeffs[1], coeffs[0]))
    plt.legend(loc='lower right')
    plt.savefig("errors.eps", bbox_inches='tight')
    plt.show()

    import ipdb; ipdb.set_trace()

if args.epsilon0:
    convergence_epsilon()
