#!/usr/bin/env python

#  TODO: Error with project (urbain, Fri 27 Apr 2018 05:23:07 PM BST)

# IMPORT MODULES {{{

import argparse
import sympy as sym
import sympy.printing as syp
# import sympy.plotting as splot
import multiprocessing
import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as las

from hermite import hermite as hm
from hermite import equations as eq
from hermite import cache as ca

from scipy.special import binom

cache = ca.cache()
sym.init_printing()

# Parse options
parser = argparse.ArgumentParser()

parser.add_argument('-c', '--config', type=str, help='Configuration file')
parser.add_argument('-b', '--beta', type=str, help='Value of β')
parser.add_argument('-t', '--theta', type=str, help='Value of θ')
parser.add_argument('-e', '--epsilon', type=str, help='Value of ε')
parser.add_argument('-m', '--mass', type=str, help='Value of m')
parser.add_argument('-p', '--plots', action='store_true',
                    help='Enable plots')
parser.add_argument('-v', '--verbose', action='store_true',
                    help='Enable verbose output')
args = parser.parse_args()

config = __import__(args.config) if args.config else __import__("config")

if args.beta:
    config.eq['β'] = sym.Rational(args.beta)
if args.theta:
    config.eq['θ'] = sym.Rational(args.theta)
if args.epsilon:
    config.eq['ε'] = sym.Rational(args.epsilon)
if args.verbose:
    config.misc['verbose'] = args.verbose
if args.plots:
    config.misc['plots'] = args.plots

config.misc['verbose'] = args.verbose

# Set library option
hm.settings.update(config.misc)

if config.misc['plots']:
    import matplotlib.pyplot as plt
    import plot

# }}}
# DATA AND PARAMETERS FOR NUMERICAL SIMULATION {{{


# Variables and function
x = eq.McKean_Vlasov.x
y = eq.McKean_Vlasov.y
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
        while value.free_symbols - {x, y, params['m'].symbol} != set():
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
    alpha = sym.sqrt(2)
    σy = Parameter(sym.symbols('σy', real=True, positive=True),
                   1/(alpha), type='equation')
    Vy = Parameter(sym.Function('Vy')(y), alpha*y*y/2,
                   type='equation')
    return {'σy': σy, 'Vy': Vy}


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
    syp.pprint(forward)
    syp.pprint(backward)

for key in params:
    params[key].value = params[key].eval()

degree = config.num['degree']
n_points_num = config.num['n_points_num']

# }}}
# DEFINE QUADRATURES {{{


def compute_quads():
    μx = params['μx'].value
    σx = params['σx'].value
    σy = params['σy'].value

    quad_num = hm.Quad.gauss_hermite(
            config.num['n_points_num'], dim=2,
            mean=[μx, 0], cov=[[σx, 0], [0, σy]])

    # Calculate limits of resolution
    band_width = np.sqrt(2) * np.sqrt(2*degree + 1)
    # x_min = μx - band_width * np.sqrt(σx)
    # x_max = μx + band_width * np.sqrt(σx)
    # y_min = 0 - band_width * np.sqrt(σy)
    # y_max = 0 + band_width * np.sqrt(σy)

    # For visualization
    nv = 200
    # bounds_x, bounds_y = band_width*np.sqrt('σx'), band_width*np.sqrt(cov_y)
    # bounds_x, bounds_y = 5*np.sqrt('σx'), 5*np.sqrt(cov_y)
    bounds_x, bounds_y = 3, 4
    quad_visu = hm.Quad.newton_cotes([nv, nv], [bounds_x, bounds_y])
    return quad_num, quad_visu


quad_num, quad_visu = compute_quads()

# }}}
# SPECTRAL METHOD FOR STATIONARY EQUATION {{{


def compute_m(series):
    quad_x = quad_num.project('x')
    series_x = series.project('x')
    x_d = quad_x.discretize('x')
    factor_d = quad_x.discretize(factor_x)
    solution_d = quad_x.eval(series_x) * factor_d
    weight_d = quad_x.discretize(quad_x.weight())
    moment0 = quad_x.integrate(solution_d/weight_d)
    moment1 = quad_x.integrate(solution_d/weight_d*x_d)
    return moment1 / moment0


def compute_moment1(m_val):
    if config.misc['verbose']:
        print("Solving the eigenvalue problem for m = " + str(m_val) + ".")
    mat_operator = r_mat + m_val * m_mat
    eig_vals, eig_vecs = cache(las.eigs)(mat_operator, k=1, which='LR')
    ground_state = np.real(eig_vecs.T[0])
    ground_state = ground_state * np.sign(ground_state[0])
    ground_series = quad_num.series(ground_state, norm=True)
    value = compute_m(ground_series)
    if config.misc['verbose']:
        print(value)
    # if config.misc['plots']:
        # fig, ax = plt.subplots(1, 1)
        # cont = quad_visu.plot(ground_series, degree, factor, ax)
        # plt.colorbar(cont, ax=ax)
        # plt.show()
    return value


m_operator = backward.diff(params['m'].symbol)
r_operator = (backward - params['m'].symbol*m_operator).cancel()
m_mat = quad_num.discretize_op(m_operator, f, degree, 2)
r_mat = quad_num.discretize_op(r_operator, f, degree, 2)

if config.misc['verbose']:
    print(hm.stats)


if args.mass:
    print(compute_moment1(float(args.mass)))
    matrix = r_mat + float(args.mass) * m_mat
    eig_vals, eig_vecs = cache(las.eigs)(matrix, k=1, which='LR')
    ground_state = np.real(eig_vecs.T[0])
    ground_state = ground_state * np.sign(ground_state[0])
    ground_series = quad_num.series(ground_state)
    ground_state_eval = quad_num.eval(ground_series)
    factor_eval = quad_num.discretize(factor)
    ground_state_eval = ground_state_eval * factor_eval
    norm = quad_num.integrate(ground_state_eval, l2=True)
    ground_state_eval = ground_state_eval / norm
    fig, ax = plt.subplots(1, 1)
    cont = quad_visu.plot(ground_series, degree, factor, ax)
    plt.colorbar(cont, ax=ax)
    plt.show()


def convergence():
    solution = eq.solve_gaussian(forward, f, [x, y])
    norm_sol = quad_num.integrate(solution, l2=True)
    assert abs(norm_sol - 1) < 1e-6
    v0 = None
    degrees = []
    errors = []
    for d in range(5, degree):
        npolys = int(binom(d + 2, d))
        sub_mat = (r_mat[0:npolys, 0:npolys]).copy(order='C')
        # sub_mat = sub_mat
        if v0 is not None:
            actual_v0 = np.zeros(npolys)
            for i in range(len(v0)):
                actual_v0[i] = v0[i]
            v0 = actual_v0
            eig_vals, eig_vecs = cache(las.eigs)(sub_mat, v0=v0, k=1, which='LR')
        else:
            eig_vals, eig_vecs = cache(las.eigs)(sub_mat, k=1, which='LR')
        ground_state = np.real(eig_vecs.T[0])
        v0 = ground_state.copy(order='C')
        ground_state = ground_state * np.sign(ground_state[0])
        ground_series = quad_num.series(ground_state)
        ground_state_eval = quad_num.eval(ground_series)
        factor_eval = quad_num.discretize(factor)
        ground_state_eval = ground_state_eval * factor_eval
        norm = quad_num.integrate(ground_state_eval, l2=True)
        ground_state_eval = ground_state_eval / norm
        solution_eval = quad_num.discretize(solution)
        error = la.norm(ground_state_eval - solution_eval, 2)
        degrees.append(d)
        errors.append(error)
    fig, ax = plt.subplots(1, 1)
    cont = quad_visu.plot(ground_series, degree, factor, ax)
    plt.colorbar(cont, ax=ax)
    plt.show()
    fig, ax = plt.subplots(1, 1)
    ax.semilogy(degrees, errors, 'k.')
    plt.show()
    # splot.contour(solution, (x, -1, 1), (y, -1, 1))


if params['Vp'].value.diff(x, x, x) == 0 and params['θ'].value == 0:
    convergence()

# asymptotic_sol_2d = sym.exp(- params['βx'].value * params['Vp'].value
#                             - params['βy'].value * params['Vy'].value)
# v0 = quad_num.transform(asymptotic_sol_2d/factor, degree, norm=True).coeffs

m_values = np.linspace(-2, 2, 41)
m_values = [0]
images = []
x_series = quad_num.transform('x', degree)

if config.misc['parallel']:
    with multiprocessing.Pool() as p:
        images = p.map(compute_moment1, m_values)
else:
    for m_val in m_values:
        images.append(compute_moment1(m_val))

# fig, ax = plt.subplots(1, 1)
# ax.plot(m_values, m_values)
# ax.plot(m_values, images)
# plt.show()
# exit(0)

# Calculate eigenvalues of largest real part
# sparse_operator = scipy.sparse.csr_matrix(mat_operator)

# }}}
# PLOTS {{{


def plot_eigenfunctions():
    fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
    axes = (ax11, ax12, ax21, ax22)
    for e_val, e_vec, ax in zip(eig_vals, eig_vecs.T, axes):
        e_vec = np.real(e_vec) * np.sign(np.real(e_vec[0]))
        series = quad_num.series(e_vec, norm=True)
        cont = quad_visu.plot(series, degree, factor, ax)
        ax.set_title("Eigenvalue: " + str(e_val))
        plt.colorbar(cont, ax=ax)
    plt.show()


def plot_hermite_functions():
    fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
    plot.plot_hf((degree, 0), quad_num, quad_visu, ax11)
    plot.plot_hf((0, degree), quad_num, quad_visu, ax12)
    plot.plot_hf((degree // 2, degree // 2), quad_num, quad_visu, ax21)
    plt.show()


def plot_ground_state():
    ground_state = np.real(eig_vecs.T[0])
    ground_state = ground_state * np.sign(ground_state[0])
    ground_series = quad_num.series(ground_state, norm=True)
    factors = {'x': factor_x, 'y': factor_y.subs(y, x)}
    plot.plot_projections(ground_series, quad_visu, factors, degree)
    plt.show()


def plot_comparison_with_asym():
    fig, ax = plt.subplots(1, 1)
    as_sol = sym.exp(- params['β'] * functions['Vp'])/factor_x
    as_series = quad_num.project('x').transform(as_sol, degree, norm=True)
    quad_visu.project('x').plot(as_series, degree, factor_x, ax)
    plt.show()


def plot_discretization_error():
    fig, ax = plt.subplots(1, 2)
    quad_visu_x = quad_visu.project('x')
    as_sol = sym.exp(- params['β'] * functions['Vp'])
    as_series = quad_num.project('x').transform(as_sol/factor_x, degree)
    norm_as = la.norm(as_series.coeffs, 2)
    as_series.coeffs = as_series.coeffs / norm_as
    quad_visu_x.plot(as_series, degree, factor_x, ax[0])
    sol_x = quad_visu_x.discretize(as_sol)
    ax[0].plot(quad_visu_x.discretize('x'), sol_x/norm_as)
    plt.show()


# plot_eigenfunctions()
# plot_hermite_functions()
plot_ground_state()
# plot_discretization_error()
# plot_comparison_with_asym()

# }}}
