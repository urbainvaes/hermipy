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
import argparse
import sympy as sym
import numpy as np
import matplotlib
import hermipy as hm
import hermipy.equations as eq

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
parser.add_argument('-l', '--lamda', type=str)
parser.add_argument('-g', '--gamma', type=float)
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

if 'DISPLAY' not in os.environ:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

# Calculate diffusion coefficients

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
α = set_param(args.alpha, r(1), 'γ') / ε**2
β = set_param(args.beta, r(1), 'β')
λ = set_param(args.lamda, r(1), 'λ') / ε

Vy = sym.cos(y)
params = {'β': β, 'α': α, 'λ': λ, 'Vy': Vy}
backward = e.backward(params)

# Numerical parameters
index_set = 'cube'
s2x, s2z, s2w = r(1, 5), r(1, 5), r(1, 5)
# s2x, s2z = r(1), r(1)
degree = args.degree if args.degree else 10
kwargs0 = {'degree': degree, 'index_set': index_set}
n_points_num = 2*degree + 1

# Calculate factor {{{
Vqx = sym.Rational(1/2)*x*x/s2x
Vqz = sym.Rational(1/2)*z*z/s2z

# Map to appropriate space
factor_y = sym.exp(- 1/2*Vy)  # Fourier
factor_x = sym.exp(- 1/2*Vqx + β/2*x*x/2)  # Hermite functions
factor_z = sym.exp(- 1/2*Vqz + β/2*z*z/2)  # Hermite functions
factor_xz = factor_x * factor_z
# }}}

# Definition of quadratures
new_q = hm.Quad.gauss_hermite
cov = [[s2x, 0], [0, s2z]]
qy = hm.Quad.fourier(2*degree + 1, dirs=[1])
qxz = new_q(2*degree + 1, dirs=[0, 2], mean=[0]*2, cov=cov, factor=factor_xz)
qf = qy * qxz  # fine
qy = hm.Quad.fourier(degree + 1, dirs=[1])
qxz = new_q(degree + 1, dirs=[0, 2], mean=[0]*2, cov=cov, factor=factor_xz)
qc = qy * qxz  # coarse

# Discretization of operator
operator = qf.discretize_op(backward, degree=degree)
rhs = qf.transform('x', degree=degree)
solution = operator.solve(rhs)
diffusion = - float(solution*rhs)
print(diffusion)
exit(0)
import ipdb; ipdb.set_trace()

# With 2 auxiliary processes
del params['α']
del params['λ']
λ1, λ2 = set_param(args.lamda, r(5), 'λ') / ε, 0
A11, A12, A21, A22 = 0, -1/ε, 1/ε, set_param(args.lamda, r(1), 'α')/ε

# Calculate factor {{{
Vqx = sym.Rational(1/2)*x*x/s2x
Vqz = sym.Rational(1/2)*z*z/s2z
Vqw = sym.Rational(1/2)*w*w/s2w

# Map to appropriate space
factor_y = sym.exp(- 1/2*Vy)  # Fourier
factor_x = sym.exp(- 1/2*Vqx + β/2*x*x/2)  # Hermite functions
factor_z = sym.exp(- 1/2*Vqz + β/2*z*z/2)  # Hermite functions
factor_w = sym.exp(- 1/2*Vqw + β/2*w*w/2)  # Hermite functions
factor = factor_x * factor_z * factor_w
# }}}
params = e.params(n_extra=2)
backward = e.backward(params)

# Fokker Planck for McKean-Vlasov equation
params.update({'Vp': Vp})
forward = equation.equation(params)

# Projections
qx, qy, qz = quad.project(0), quad.project(1), quad.project(2)
wx = qx.factor * qx.factor / qx.position.weight()
wy = qy.factor * qy.factor / qy.position.weight()
wz = qz.factor * qz.factor / qz.position.weight()
qxy, qxz, qyz = quad.project([0, 1]), quad.project([0, 2]), quad.project([1, 2])

# Integral operators
Ix = qx.transform(wx, degree=degree, index_set=index_set)
Iy = qy.transform(wy, degree=degree, index_set=index_set)
Iz = qz.transform(wz, degree=degree, index_set=index_set)

# Moment operators
mx1 = qx.transform(wx * x, degree=degree, index_set=index_set)
my1 = qy.transform(wy * y, degree=degree, index_set=index_set)
mz1 = qz.transform(wz * z, degree=degree, index_set=index_set)

if args.interactive:
    plt.ion()
    fig, ax = plt.subplots(2, 2)


def plot(t, βplot=β):

    # Retrieve degree
    d = t.degree

    # Integral operators
    Ix_d = Ix.subdegree(d)
    Iy_d = Iy.subdegree(d)
    Iz_d = Iz.subdegree(d)

    # Projection on x-q subspace
    txy, txz = Iz_d*t, Iy_d*t

    # Projections
    tx, ty, tz = Iy_d*(Iz_d*t), Ix_d*(Iz_d*t), Ix_d*(Iy_d*t)

    # Moments
    mx = str((mx1.subdegree(d)*tx).coeffs[0])
    my = str((my1.subdegree(d)*ty).coeffs[0])
    mz = str((mz1.subdegree(d)*tz).coeffs[0])

    ax[0][0].clear()
    plot_position = True
    field = txz if plot_position else txy
    qxz.plot(field, ax=ax[0][0], bounds=False, vmin=0, extend='min')

    ax[0][1].clear()
    qx.plot(tx, ax=ax[0][1], title="1st moment - X: " + mx)
    if not args.bifurcation:
        qx.plot(ux, ax=ax[0][1])

    ax[1][0].clear()
    qy.plot(uy, ax=ax[1][0])
    qy.plot(ty, ax=ax[1][0], title="1st moment - P: " + my)

    ax[1][1].clear()
    qz.plot(uz, ax=ax[1][1])
    qz.plot(tz, ax=ax[1][1], title="1st moment - Q: " + mz)

    plt.draw()
    plt.pause(0.1)
