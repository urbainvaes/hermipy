#!/usr/bin/env python

# IMPORT MODULES {{{

import hashlib
import sympy as sym
import sympy.printing as syp
import sympy.plotting as splot
import numpy as np
import numpy.linalg as la
# import scipy.sparse
import scipy.sparse.linalg as las
import matplotlib.pyplot as plt

import plot
import equation
import config
from libhermite import hermite as hm

sym.init_printing()

# }}}
# DATA AND PARAMETERS FOR NUMERICAL SIMULATION {{{

# Set library option
hm.rc['cache'] = config.glob['cache']

# Variables and function
x, y, f = equation.x, equation.y, equation.f

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
    for key, symbol in equation.forward_params().items():
        value = config.eq[key] if key in config.eq else symbol
        eq_params[key] = Parameter(symbol, value, type='equation')
    return eq_params


def vy():
    alpha = sym.sqrt(2/params['βy'].symbol)
    σy = Parameter(sym.symbols('σy', real=True, positive=True),
                   1/(alpha*params['βy'].symbol),
                   type='equation')
    Vy = Parameter(sym.Function('Vy')(y),
                   y*y/(2*params['βy'].symbol*σy.symbol),
                   type='equation')
    return {'σy': σy, 'Vy': Vy}


def vq():
    r = sym.Rational
    symb = sym.symbols
    μx = Parameter(symb('μx', real=True), config.num['μx'])
    σx = Parameter(symb('σx', real=True, positive=True), config.num['σx'])
    Vq = Parameter(sym.Function('Vq')(x),
                   r(1, 2)*(x - μx.symbol)*(x - μx.symbol)/σx.symbol)
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
    return equation.forward(eq_params)


def factors(symbolic, λ):
    βx = params['βx'].value if symbolic == 0 else params['βx'].symbol
    βy = params['βy'].value if symbolic == 0 else params['βy'].symbol
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
    factor_x = sym.exp(-βx*(λ*Vq + (1-λ)*Vp))
    factor_y = sym.exp(-βy*Vy)
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
forward = forward_op(config.glob['symbolic'])

syp.pprint(forward)

# Mapped backward operator
factor_x, factor_y, factor = factors(config.glob['symbolic'], config.num['λ'])
backward = equation.map_operator(forward, factor)

forward, backward = evaluate([forward, backward])
factor_x, factor_y, factor = evaluate([factor_x, factor_y, factor])

syp.pprint(backward)

for key in params:
    params[key].value = params[key].eval()

degree = config.num['degree']
n_points_num = config.num['n_points_num']

# print(r_operator)

# }}}
# DEFINE QUADRATURES {{{

if params['Vp'].value.diff(x, x, x) == 0:
    solution = equation.solve_gaussian(forward)
    splot.plot3d(solution, (x, -1, 1), (y, -1, 1))

m_operator = backward.diff(params['m'].symbol)
r_operator = (backward - params['m'].symbol*m_operator).cancel()


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


m_mat = quad_num.discretize_op(m_operator, f, degree, 2)
r_mat = quad_num.discretize_op(r_operator, f, degree, 2)


# asymptotic_sol_2d = sym.exp(- params['βx'].value * params['Vp'].value
#                             - params['βy'].value * params['Vy'].value)
# v0 = quad_num.transform(asymptotic_sol_2d/factor, degree, norm=True).coeffs

m_values = np.linspace(-1, 1, 11)
m_values = [0]
images = []
x_series = quad_num.transform('x', degree)
print(x_series)

for m_val in m_values:
    mat_operator = r_mat + m_val * m_mat
    print("Solving the eigenvalue problem...")
    eig_vals, eig_vecs = hm.cache(las.eigs)(mat_operator, k=1, which='LR')
    ground_state = np.real(eig_vecs.T[0])
    ground_series = quad_num.series(ground_state, norm=True)
    value = compute_m(ground_series)
    images.append(value)
    print(value)

fig, ax = plt.subplots(1, 1)
ax.plot(m_values, m_values)
ax.plot(m_values, images)

# Calculate eigenvalues of largest real part
# sparse_operator = scipy.sparse.csr_matrix(mat_operator)

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
    ground_series = quad_num.series(ground_state, norm=True)
    factors = {'x': factor_x, 'y': factor_y.subs(y, x)}
    plot.plot_projections(ground_series, quad_visu, factors, degree)
    plt.show()


def plot_comparison_with_asym():
    fig, ax = plt.subplots(1, 1)
    as_sol = sym.exp(- params['βx'] * functions['Vp'])/factor_x
    as_series = quad_num.project('x').transform(as_sol, degree, norm=True)
    quad_visu.project('x').plot(as_series, degree, factor_x, ax)
    plt.show()


def plot_discretization_error():
    fig, ax = plt.subplots(1, 2)
    quad_visu_x = quad_visu.project('x')
    as_sol = sym.exp(- params['βx'] * functions['Vp'])
    as_series = quad_num.project('x').transform(as_sol/factor_x, degree)
    norm_as = la.norm(as_series.coeffs, 2)
    as_series.coeffs = as_series.coeffs / norm_as
    quad_visu_x.plot(as_series, degree, factor_x, ax[0])
    sol_x = quad_visu_x.discretize(as_sol)
    ax[0].plot(quad_visu_x.discretize('x'), sol_x/norm_as)
    plt.show()


plot_eigenfunctions()
# plot_hermite_functions()
plot_ground_state()
# plot_discretization_error()
# plot_comparison_with_asym()

# }}}
