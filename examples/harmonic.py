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
import numpy.linalg as la
import matplotlib
import hermipy as hm
import hermipy.equations as eq
import hermipy.core as core

hm.settings['tensorize'] = True
hm.settings['sparse'] = True
hm.settings['trails'] = False
hm.settings['cache'] = True

# Process arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--interactive', action='store_true')
parser.add_argument('-tcd', '--convergence_degree', action='store_true')
parser.add_argument('-tce', '--convergence_epsilon', action='store_true')
parser.add_argument('-e', '--epsilon', type=float)
parser.add_argument('-b', '--beta', type=float)
parser.add_argument('-t', '--theta', type=float)
parser.add_argument('-g', '--gamma', type=float)
parser.add_argument('-m', '--mass', type=float)
args = parser.parse_args()

# Matplotlib configuration
matplotlib.rc('font', size=14)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)

if 'DISPLAY' not in os.environ:
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

# Dimension of the problem
dim = 3

# Equation parameters
index_set = 'rectangle'
r, s = sym.Rational, sym.symbols
equation = eq.McKean_Vlasov_harmonic_noise
x, y, z, f = equation.x, equation.y, equation.z, equation.f
β, ε, γ, θ, m = r(5), r(1, 5), 0, 0, 0

# For convergence_epsilon
ε, γ = s('ε'), s('γ')

params = {'β': β, 'ε': ε, 'γ': γ, 'θ': θ, 'm': m}
Vp, β = x**4/4 - x**2/2, params['β']

# Numerical parameters
s2x, s2y, s2z, degree = r(1, 20), r(1, 5), r(1, 5), 50
n_points_num = 2*degree + 1

# Potential for approximation
Vqx = sym.Rational(1/2)*x*x/(β*s2x)
Vqy = sym.Rational(1/2)*y*y/s2y
Vqz = sym.Rational(1/2)*z*z/s2z

# Fokker Planck for McKean-Vlasov equation
params.update({'Vp': Vp})
forward = equation.equation(params)

# Map to appropriate space
factor_x = sym.exp(- β / 2 * (Vqx + Vp))
factor_pq = sym.exp(- 1/2 * (y*y/2 + Vqy + z*z/2 + Vqz))
factor = factor_x * factor_pq

# Calculation of the solution
new_q = hm.Quad.gauss_hermite
cov = [[s2x, 0, 0], [0, s2y, 0], [0, 0, s2z]]
quad = new_q(n_points_num, dim=3, mean=[0]*3, cov=cov, factor=factor)
min_quad = new_q(n_points=degree+1, dim=3, mean=[0]*3, cov=cov, factor=factor)

# Projections
qx, qy, qz = quad.project(0), quad.project(1), quad.project(2)
wx = qx.factor * qx.factor / qx.position.weight()
wy = qy.factor * qy.factor / qy.position.weight()
wz = qz.factor * qz.factor / qz.position.weight()
qxy, qxz = quad.project([0, 1]), quad.project([0, 2])

# Integral operators
Ix = qx.transform(wx, degree=degree, index_set=index_set)
Iy = qy.transform(wy, degree=degree, index_set=index_set)
Iz = qz.transform(wz, degree=degree, index_set=index_set)

# Moment operators
mx1 = qx.transform(wx * x, degree=degree, index_set=index_set)
my1 = qy.transform(wy * y, degree=degree, index_set=index_set)
mz1 = qz.transform(wz * z, degree=degree, index_set=index_set)

# Marginals in white noise case
ux = sym.exp(-β*Vp) / qx.integrate(sym.exp(-β*Vp), flat=True)
uy = sym.exp(-y*y/2) / qy.integrate(sym.exp(-y*y/2), flat=True)
uz = sym.exp(-z*z/2) / qz.integrate(sym.exp(-z*z/2), flat=True)

if args.interactive:
    plt.ion()
    fig, ax = plt.subplots(2, 2)


def plot(t):

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
    qxz.plot(txz, ax=ax[0][0], bounds=False, vmin=0, extend='min')
    qxy.plot(txy, ax=ax[0][0], bounds=False, vmin=0, extend='min')

    ax[0][1].clear()
    qx.plot(ux, ax=ax[0][1])
    qx.plot(tx, ax=ax[0][1], title="1st moment - X: " + mx)

    ax[1][0].clear()
    qy.plot(uy, ax=ax[1][0])
    qy.plot(ty, ax=ax[1][0], title="1st moment - P: " + my)

    ax[1][1].clear()
    qz.plot(uz, ax=ax[1][1])
    qz.plot(tz, ax=ax[1][1], title="1st moment - Q: " + mz)

    plt.draw()
    plt.pause(0.1)


def convergence_degree():

    eye = quad.varf('1', degree=degree, index_set=index_set)
    dt, Ns, scheme = 2**-9, int(1e4), "backward"
    dmin, d, degrees = 10, degree, []
    errors, mins = [], []

    # Initial condition
    t = quad.transform(ux*uy*uz, degree=degree, index_set=index_set)

    while d >= dmin:
        r_mat, eye = mat.subdegree(d), eye.subdegree(d)
        t = t.subdegree(d)

        # Integral operators
        Ix_d = Ix.subdegree(d)
        Iy_d = Iy.subdegree(d)
        Iz_d = Iz.subdegree(d)

        Δ = np.inf
        for i in range(Ns):

            if args.interactive:
                plot(t)

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
            integral = (Ix_d*(Iy_d*(Iz_d*new_t))).coeffs[0]
            new_t = new_t / integral

            # Error
            Δ = ((t - new_t)*(t - new_t)).coeffs[0]/dt

            # Time adaptation
            threshold, dt_max = .01, 64
            if Δ*dt > 2*threshold:
                dt = dt / 2.
                continue
            elif Δ*dt < threshold and dt < dt_max:
                dt = dt * 2.
            t = new_t

            if Δ < 1e-15:

                if d == degree:
                    t_exact = t

                else:
                    error_series = t.subdegree(degree) - t_exact
                    error_series_x = Iy*(Iz*error_series)
                    error = min_quad.norm(error_series, n=1, flat=True)
                    error_x = qx.norm(error_series_x, n=1, flat=True)
                    min = np.min(min_quad.eval(error_series))
                    mins.append(abs(min))
                    errors.append(error)
                    degrees.append(d)
                    print(min, error, error_x)

                if args.interactive:
                    plot(t)

                d = d - 1
                break

    cond = np.asarray(degrees)*0 + 1
    xplot, yplot = np.extract(cond, degrees), np.extract(cond, errors)
    fig, ax = plt.subplots()
    ax.semilogy(xplot, yplot, 'b.', label="$\\|\\rho^{{50}} - \\rho^d\\|_1$")
    coeffs = np.polyfit(xplot, np.log10(yplot), 1)
    ax.semilogy(xplot, 10**coeffs[1] * 10**(coeffs[0]*xplot), 'b-')
    yplot = np.extract(cond, mins)
    ax.semilogy(xplot, yplot, 'r.', label="$|\\min \\, \\rho^d(x, p, q)|$")
    coeffs = np.polyfit(xplot, np.log10(yplot), 1)
    ax.semilogy(xplot, 10**coeffs[1] * 10**(coeffs[0]*xplot), 'r-')
    ax.set_xlabel("$d$")
    plt.legend(loc='upper right')
    plt.savefig("errors.eps", bbox_inches='tight')


if args.convergence_degree:
    convergence_degree()


def convergence_epsilon():

    # Initial condition
    t = quad.transform(ux*uy*uz, degree=degree, index_set=index_set)

    eye = quad.varf('1', degree=degree, index_set=index_set)
    dt, Ns, scheme = 2**-9, int(1e4), "backward"

    # Calculate white noise solution
    kwargs = {'degree': degree, 'index_set': index_set}
    forward0x = eq.Fokker_Planck_1d.equation({'Vp': Vp, 'β': β})
    var0x = qx.discretize_op(forward0x, **kwargs)
    _, [tx0] = var0x.eigs(k=1, which='LR')
    forward0 = forward.subs(γ, 1).subs(ε, 1)
    var0 = quad.discretize_op(forward0, **kwargs)
    _, [t30] = var0.eigs(k=1, which='LR')
    tx0, t30 = tx0 / float(Ix*tx0), t30 / float(Ix*(Iy*(Iz*t30)))

    εs = [2**(-i/4) for i in range(30)]
    e3, ex, t3, tx = [], [], [], []

    for iε, εi in enumerate(εs):

        forward_ε = forward.subs(γ, 0).subs(ε, εi)
        r_mat = quad.discretize_op(forward_ε, **kwargs)

        Δ, Δx = np.inf, np.inf
        for i in range(Ns):

            if args.interactive:
                plot(t)

            print("ε: " + str(εi) + ", i: " + str(i) + ", dt: " + str(dt)
                  + ", Δ: " + str(Δ) + ", Δx: " + str(Δx))

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
            new_t = new_t / float(Ix*(Iy*(Iz*new_t)))

            # Error
            Δ = ((t - new_t)*(t - new_t)).coeffs[0]/dt

            # Time adaptation
            threshold, dt_max = .01, 64
            if Δ*dt > 2*threshold:
                dt = dt / 2.
                continue
            elif Δ*dt < threshold and dt < dt_max:
                dt = dt * 2.
            t = new_t

            if Δ < 1e-16:

                if args.interactive:
                    plot(t)

                t3.append(t)
                tx.append(Iy*(Iz*t))
                error3, errorx = t30 - t3[-1], tx0 - tx[-1]
                e3.append(min_quad.norm(error3, n=1, flat=True))
                ex.append(qx.norm(errorx, n=1, flat=True))

                if iε > 1:
                    logε = np.log2(εs[0:iε+1])
                    coeffs = np.polyfit(logε, np.log2(e3), 1)
                    print(coeffs)
                    coeffs_x = np.polyfit(logε, np.log2(ex), 1)
                    print(coeffs_x)

                break

    t30, tx0 = t3[-1], tx[-1]
    ex = [qx.norm(tx0 - txi, n=1, flat=True) for txi in tx]
    e3 = [min_quad.norm(t30 - t3i, n=1, flat=True) for t3i in t3]
    εs, ex, e3 = εs[0:-1], ex[0:-1], e3[0:-1]

    fig, ax = plt.subplots()
    cmap = matplotlib.cm.get_cmap('viridis_r')
    for i, (εi, txε) in enumerate(zip(εs, tx)):
        if i % 2 is not 0:
            continue
        label = "$\\varepsilon = 2^{{ -{} }} $".format(i/4)
        kwargs = {'color': cmap(εi), 'label': label}
        qx.plot(txε, ax=ax, **kwargs)
    ax.set_title("")
    plt.legend()
    plt.savefig("convergence-epsilon.eps", bbox_inches='tight')

    fig, ax = plt.subplots()
    xplot, yplot1 = logε = np.asarray(εs), np.asarray(e3)
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


if args.convergence_epsilon:
    convergence_epsilon()

# Discretization of the operator
degrees = list(range(5, degree + 1))
mat = quad.discretize_op(forward, degrees[-1], index_set=index_set)
solutions, v0, eig_vec = [], None, None

if args.interactive:
    plt.ion()
    fig, ax = plt.subplots(2, 2)

for d in degrees:
    print(d)
    npolys = core.iterator_size(dim, d, index_set=index_set)
    if d is not degrees[0]:
        v0 = np.zeros(npolys)
        for i in range(len(eig_vec.coeffs)):
            v0[i] = eig_vec.coeffs[i]
    sub_mat = mat.subdegree(d)
    eig_vals, [eig_vec] = sub_mat.eigs(k=1, v0=v0, which='LR')
    t = eig_vec * np.sign(eig_vec.coeffs[0])

    t = t / quad.integrate(t, flat=True, tensorize=False)
    solutions.append(t)

    if args.interactive:

        Ix_d = Ix.subdegree(d)
        Iy_d = Iy.subdegree(d)
        Iz_d = Iz.subdegree(d)

        # Projection on x-q subspace
        integral = (Ix_d*(Iy_d*(Iz_d*t))).coeffs[0]
        txy, txz = Iz_d*t, Iy_d*t

        # Projections
        tx, ty, tz = Iy_d*(Iz_d*t), Ix_d*(Iz_d*t), Ix_d*(Iy_d*t)

        # Moments
        mx = str((mx1.subdegree(d)*tx).coeffs[0])
        my = str((my1.subdegree(d)*ty).coeffs[0])
        mz = str((mz1.subdegree(d)*tz).coeffs[0])

        ax[0][0].clear()
        qxz.plot(txz, ax=ax[0][0], bounds=False, vmin=0, extend='min')
        qxy.plot(txy, ax=ax[0][0], bounds=False, vmin=0, extend='min')

        ax[0][1].clear()
        qx.plot(ux, ax=ax[0][1])
        qx.plot(tx, ax=ax[0][1], title="1st moment - X: " + mx)

        ax[1][0].clear()
        qy.plot(uy, ax=ax[1][0])
        qy.plot(ty, ax=ax[1][0], title="1st moment - P: " + my)

        ax[1][1].clear()
        qz.plot(uz, ax=ax[1][1])
        qz.plot(tz, ax=ax[1][1], title="1st moment - Q: " + mz)

        plt.draw()
        plt.pause(.01)


finest = ground_state
finest_eval = solutions[-1]

# Plot of the finest solution
fig, ax = plt.subplots(1, 1)
quad.plot(finest, factor, ax=ax)
plt.show()

# Associated errors
errors, degrees = [], degrees[0:-1]
for sol in solutions[0:-1]:
    error = quad.norm(sol - finest_eval, flat=True)
    errors.append(error)
    print(error)

log_errors = np.log(errors)
poly_approx = np.polyfit(degrees, log_errors, 1)
errors_approx = np.exp(np.polyval(poly_approx, degrees))
error = la.norm(log_errors - np.log(errors_approx), 2)

plt.semilogy(degrees, errors, 'k.')
plt.semilogy(degrees, errors_approx)
plt.show()
