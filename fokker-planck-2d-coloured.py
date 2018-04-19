#!/usr/bin/env python

# IMPORT MODULES {{{

import math
import sympy as sym
import sympy.printing as syp
import sympy.plotting as splot
import numpy as np
import numpy.linalg as la
import scipy.sparse
import scipy.sparse.linalg as las
import matplotlib.pyplot as plt

import plot
import equation
from libhermite import hermite as hm
from config import glob, params, sym_params, functions, numerics

sym.init_printing()

# }}}
# DATA AND PARAMETERS FOR NUMERICAL SIMULATION {{{

config_dict = {**params, **functions, **numerics}
hash_config = hash(frozenset(config_dict.items()))

# Short-hand notation
r = sym.Rational

# Variables and function
x, y, f = equation.x, equation.y, equation.f

# Parameters of approximating potential
def gaussian(mean, var):
    return r(.5)*(x - mean)*(x - mean)/(var)


sym_params['μx'] = sym.symbols('μx', real=True)
sym_params['σx'] = sym.symbols('σx', real=True, positive=True)
functions['Vq'] = gaussian(sym_params['μx'], sym_params['σx'])/sym_params['βx']

# Real and functional parameters
sym_functions = {}
equation.fill_params(sym_params, sym_functions)
sp, p = sym_params, params

# Parameters of potential in y equation
alpha = sym.sqrt(2/sp['βy'])
cov_y = 1/(alpha*sp['βy'])
functions['Vy'] = y*y/(2*sp['βy']*cov_y)

degree = numerics['degree']
n_points_num = numerics['n_points_num']

# Use cache
cache = glob['cache']
fastsym = True

# }}}
# ABSTRACT SYMBOLIC CALCULATIONS {{{

# Fill missing params with symbolic values
equation.fill_params(params, functions)

# Conversion lists
all_params = { **params, **numerics }
real_params = [(sym_params[k], all_params[k]) for k in sym_params]
functional_params = [(sym_functions[k], functions[k]) for k in functions]

# Evaluate parameters in functions
for key in functions:
    functions[key] = functions[key].subs(real_params)
    functions[key] = functions[key].subs(real_params)

# Operator
param_factor = r(.5)
forward, backward, factor_x, factor_y, factor = \
        equation.operators(params, functions, param_factor)

# Print operator to output
# syp.pprint(backward)

cov_y = float(cov_y.subs(real_params))

for key in params:
    if key is not 'm':
        params[key] = float(params[key])

for key in numerics:
    numerics[key] = float(numerics[key])

if functions['Vp'].diff(x, x, x) == 0:
    solution = equation.solve_gaussian(forward)
    splot.plot3d(solution.subs(real_params), (x, -1, 1), (y, -1, 1))

m_operator = backward.diff(sp['m'])
r_operator = (backward - sp['m']*m_operator).cancel()

# }}}
# DEFINE QUADRATURE FOR NUMERICS {{{

# For numerics
new_q = hm.Quad.gauss_hermite
mean_xy, cov_xy = [numerics['μx'], 0], [[numerics['σx'], 0], [0, cov_y]]
quad_num = new_q(n_points_num, dim=2, mean=mean_xy, cov=cov_xy)

# Calculate limits of resolution
band_width = np.sqrt(2) * np.sqrt(2*degree + 1)
x_min = numerics['μx'] - band_width * np.sqrt(numerics['σx'])
x_max = numerics['μx'] + band_width * np.sqrt(numerics['σx'])
y_min = 0 - band_width * np.sqrt(cov_y)
y_max = 0 + band_width * np.sqrt(cov_y)

# }}}
# DEFINE QUADRATURES FOR VISUALIZATION {{{

nv = 200
# bounds_x, bounds_y = band_width*np.sqrt('σx'), band_width*np.sqrt(cov_y)
# bounds_x, bounds_y = 5*np.sqrt('σx'), 5*np.sqrt(cov_y)
bounds_x, bounds_y = 3, 4

quad_visu = hm.Quad.newton_cotes([nv, nv], [bounds_x, bounds_y])

# }}}
# CALCULATE ASYMPTOTIC SOLUTION {{{
# }}}
# SPECTRAL METHOD FOR STATIONARY EQUATION {{{

def compute_with_cache(cache):

    args_m = [m_operator, f, degree, 2]
    args_r = [r_operator, f, degree, 2]

    str_args = [str(e).encode('utf-8') for e in args_m + args_r]
    hash_args = '-'.join(hashlib.md5(s).hexdigest() for s in str_args)
    hash_config = hashlib.md5(hash_args.encode('utf-8')).hexdigest()

    def call_discretize():
        m_mat =  quad_num.discretize_op(*args_m)
        r_mat =  quad_num.discretize_op(*args_r)
        return m_mat, r_mat

    try:
        m_mat_cache = np.load('cache/m_mat-' + str(hash_config) + '.npy')
        r_mat_cache = np.load('cache/r_mat-' + str(hash_config) + '.npy')

    except IOError:
        m_mat, r_mat = call_discretize()
        np.save('cache/m_mat-' + str(hash_config), m_mat)
        np.save('cache/r_mat-' + str(hash_config), r_mat)
        return m_mat, r_mat

    if not cache:
        m_mat, r_mat = call_discretize()
        assert la.norm(m_mat - m_mat_cache, 2) < 1e-10
        assert la.norm(r_mat - r_mat_cache, 2) < 1e-10

    return m_mat, r_mat

m_mat, r_mat = compute_with_cache(glob['cache'])

# m_values = np.linspace(-1, 1, 11)


m_num = 0
# for i in range(10):
#     total_mat = r_mat + m_num * m_mat


# Calculate eigenvalues of largest real part
print("Solving the eigenvalue problem...")
mat_operator = r_mat
sparse_operator = scipy.sparse.csr_matrix(mat_operator)
asymptotic_sol_2d = sym.exp(- params['βx'] * functions['Vp']
                            - params['βy'] * functions['Vy'])
v0 = quad_num.transform(asymptotic_sol_2d/factor, degree, norm=True).coeffs
eigen_values, eigen_vectors = las.eigs(mat_operator, v0=v0, k=4, which='LR', ncv=10)

## PLOTS {{{

def plot_eigenfunctions():
    fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
    axes = (ax11, ax12, ax21, ax22)
    for e_val, e_vec, ax in zip(eigen_values, eigen_vectors.T, axes):
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
    ground_state = np.real(eigen_vectors.T[0])
    ground_series = quad_num.series(ground_state, norm=True)
    plot.plot_projections(ground_series, quad_visu, factor, degree)
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
# plot_ground_state()
# plot_discretization_error()
# plot_comparison_with_asym()

# }}}
