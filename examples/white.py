#1 Copyright (C) 2018 Urbain Vaes

import os
import argparse
import sympy as sym
import numpy as np
import matplotlib
import hermipy as hm
import hermipy.equations as eq
import matplotlib.pyplot as plt
import scipy.io

hm.settings['tensorize'] = True
hm.settings['sparse'] = True
hm.settings['trails'] = False
hm.settings['cache'] = True

# Process arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', '--directory', type=str)
parser.add_argument('-i', '--interactive', action='store_false')
parser.add_argument('-b', '--beta', type=str)
parser.add_argument('-m0', '--mass0', type=str)
parser.add_argument('-m', '--mass', type=str)
parser.add_argument('-t', '--theta', type=str)
parser.add_argument('-d', '--degree', type=int)
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

# Dimension of the problem
dim = 1

# Equation parameters
r, s = sym.Rational, sym.symbols
equation = eq.McKean_Vlasov_white
x, f = equation.x, equation.f


def set_param(arg, default, symbol):
    if arg and r(arg) < 0:
        return s(symbol, real=True)
    if arg:
        return r(arg)
    return default


# Fokker Planck equation
β = set_param(args.beta, r(1), 'β')
m0 = set_param(args.mass0, r(0), 'm0')
θ = set_param(args.theta, r(1, 1), 'θ')
m, Vp = s('m', real=True), x**4/4 - x**2/2
params = {'β': β, 'm': m, 'θ': θ, 'Vp': Vp}
forward = equation.equation(params)

# Load data
data_file = 'bdata/sol_beta_{}_mean_{}.mat'.format(str(β), str(float(m0)))
try:
    data = scipy.io.loadmat(data_file)
    time, rho, xs = data['t'].T[0][1:], data['p'], data['x'][0]
except IOError:
    exit(0)

# Define quadrature for comparison
quad_comparison = hm.Quad(nodes=[xs], weights=[xs*0], types=["visu"])

sx, mx = r(1, 5), r(0)
degree = args.degree if args.degree else 100
kwargs0 = {'degree': degree}
n_points_num = 2*degree + 1

# Potential for approximation
Vq = sym.Rational(1/2)*x*x/sx
factor = sym.exp(-Vq/2)

# Calculation of the solution
new_q = hm.Quad.gauss_hermite
mean, cov = [mx], [[sx]]
quad = new_q(n_points_num, mean=mean, cov=cov, factor=factor)

# Integral and moment operator
w = quad.factor * quad.factor / quad.position.weight()
In = quad.transform(w, degree=degree)
m1 = quad.transform(w * x, degree=degree)

# Degrees
m_operator = forward.diff(m)
r_operator = (forward - m*m_operator).cancel()
r_mat = quad.discretize_op(r_operator, **kwargs0)
m_mat = quad.discretize_op(m_operator, **kwargs0)

# Initial condition
u = sym.exp(-(x-m0)*(x-m0)/2)
u = u / quad.integrate(u, flat=True)
t = quad.transform(u, degree=degree)
m = float(m1*t)

if args.interactive:
    plt.ion()
    fig, ax = plt.subplots(2)


def plot(i, t, title=None):

    if not args.interactive:
        return

    ax[0].clear()
    quad_comparison.plot(t, ax=ax[0], bounds=False, title=title)
    quad_comparison.plot(rho[i], ax=ax[0], bounds=False, title="Solutions")

    ax[1].clear()
    evalualion = quad_comparison.eval(t)
    quad_comparison.plot(rho[i] - evalualion, ax=ax[1], bounds=False, title="Error")
    # t.plot(ax=ax[1])

    plt.draw()
    plt.pause(0.01)


eye = quad.varf('1', degree=degree)
steps = np.diff(time)
scheme = 'crank'

for i, dt in enumerate(steps):

    # Error
    Δ = quad_comparison.eval(t) - rho[i]
    Δ = np.sum(np.abs(Δ))*(xs[-1] - xs[0])/(len(Δ) - 1)

    if args.interactive and i % 10 == 0:
        plot(i, t)

    print("i: {}, t: {}, m: {}, Δ: {}".format(i, time[i], m, Δ))

    mat = r_mat + m*m_mat
    if scheme == 'backward':
        total_op = eye - dt*mat
        new_t = total_op.solve(t)

    if scheme == 'crank':
        crank_left = eye - dt*mat*(1/2)
        crank_right = eye + dt*mat*(1/2)
        new_t = crank_left.solve(crank_right(t))

    # Normalization
    m = float(m1*new_t)
    t = new_t / float(In*new_t)
