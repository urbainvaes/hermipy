# Copyright (C) 2018 Urbain Vaes

import os
import argparse
import sympy as sym
import numpy as np
import hermipy as hm
import hermipy.equations as eq
import scipy.io
import scipy.integrate
import matplotlib

hm.settings['tensorize'] = True
hm.settings['sparse'] = True
hm.settings['trails'] = False
hm.settings['cache'] = False

# Process arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', '--directory', type=str)
parser.add_argument('-b', '--beta', type=str)
parser.add_argument('-m0', '--mass0', type=str)
parser.add_argument('-m', '--mass', type=str)
parser.add_argument('-t', '--theta', type=str)
parser.add_argument('-d', '--degree', type=int)
parser.add_argument('-i', '--interactive', action='store_true')
parser.add_argument('-ta', '--test_animate', action='store_true')
parser.add_argument('-tac', '--test_animate_coefficients', action='store_true')
parser.add_argument('-tp', '--test_plots', action='store_true')
parser.add_argument('-te', '--test_eigs', action='store_true')
parser.add_argument('-tc', '--test_convergence', action='store_true')
parser.add_argument('--method', type=str)
args = parser.parse_args()

if not args.interactive:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Directory for output files
dir = ""
if args.directory:
    dir = args.directory + "/"
    os.makedirs(dir, exist_ok=True)

# Matplotlib configuration
matplotlib.rc('font', size=20)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)
matplotlib.rc('figure', figsize=(14, 8))

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

sx, mx = r(1, 10), r(0)
degree = args.degree if args.degree else 100
kwargs0 = {'degree': degree}
n_points_num = 2*degree + 1

# Potential for approximation
Vq = sym.Rational(1/2)*x*x/sx
# factor = sym.exp(-Vq/2 - β*Vp/2)
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

if args.test_eigs:
    degrees, eigvals = list(range(4, degree)), []
    for d in degrees:
        e, v = r_mat.subdegree(d).eigs(k=1, which='LR')
        eigvals.append(e)
        print(e)
    fig, ax = plt.subplots()
    ax.plot(degrees, eigvals, 'k.')
    ax.set_yscale('symlog', linthreshy=1e-5)
    ax.set_xlabel('Number of Hermite functions')
    ax.set_title('$\\lambda_0$')
    plt.savefig('eigen.pdf', bbox_inches='tight')
    plt.show()

# Load data
data_file = 'bdata/sol_beta_{}_mean_{}.mat'.format(str(β), str(float(m0)))
try:
    data = scipy.io.loadmat(data_file)
    time, rho, xs = data['t'].T[0][1:], data['p'], data['x'][0]
except IOError:
    exit(0)

subsample = 1
time = [time[i] for i in range(0, len(time) - 1, subsample)]
rho = [rho[i] for i in range(0, len(rho) - 1, subsample)]

# Define quadrature for comparison
quad_comparison = hm.Quad(nodes=[xs], weights=[xs*0], types=["visu"])


def error_Linf(result):
    errors = np.zeros(len(result))
    for i, t in enumerate(result):
        Δ = quad_comparison.eval(t) - rho[i]
        errors[i] = np.sum(np.abs(Δ))*(xs[-1] - xs[0])/(len(Δ) - 1)
    return np.max(errors)


def solve(method, subdegree=degree):

    # Subdegree
    sub_r_mat = r_mat.subdegree(subdegree)
    sub_m_mat = m_mat.subdegree(subdegree)
    sub_In = In.subdegree(subdegree)
    sub_m1 = m1.subdegree(subdegree)

    # Initial condition
    u = sym.exp(-(x-m0)*(x-m0)/2)
    u = u / quad.integrate(u, flat=True)
    t = quad.transform(u, degree=subdegree)
    m = float(sub_m1*t)

    if method == "ode45":

        def dfdt(t, y):
            m = np.dot(sub_m1.coeffs, y) / np.dot(sub_In.coeffs, y)
            mat = sub_r_mat.matrix + m*sub_m_mat.matrix
            return_vector = mat.dot(y)
            return return_vector
        result = scipy.integrate.solve_ivp(dfdt, [0, 5], t.coeffs, 'RK45',
                                           t_eval=time, max_step=.01,
                                           atol=1e-9, rtol=1e-9)
        result = [quad.series(y) for y in result.y.T]

    elif method == "semi_explicit":

        eye = quad.varf('1', degree=subdegree)
        steps = np.diff(time)
        scheme = 'crank'
        result = [t]

        for i, dt in enumerate(steps):

            print("i: {}, t: {}, m: {}".format(i, time[i], m))

            mat = sub_r_mat + m*sub_m_mat
            if scheme == 'backward':
                total_op = eye - dt*mat
                new_t = total_op.solve(t)

            if scheme == 'crank':
                crank_left = eye - dt*mat*(1/2)
                crank_right = eye + dt*mat*(1/2)
                new_t = crank_left.solve(crank_right(t))

            # Normalization
            m = float(sub_m1*new_t) / float(sub_In*new_t)

            # t = new_t / float(In*new_t)
            t = new_t
            result.append(t)

    return result


if args.test_convergence:

    try:
        degrees = np.load("data/comparison_carillo/comparison_carrillo_degrees_ode45_80.npy")
        errors_ode45 = np.load("data/comparison_carillo/comparison_carrillo_errors_ode45_80.npy")
        errors_semi_explicit = np.load("data/comparison_carillo/comparison_carrillo_errors_semi_explicit_80.npy")

        fig, ax = plt.subplots()
        plt.semilogy(degrees, errors_ode45, '.', label='RK45')
        plt.semilogy(degrees, errors_semi_explicit, '.', label='Semi-implicit')
        ax.set_xlabel("Number of Hermite functions")
        ax.set_title("Difference (in the $L^\\infty(0, T; L^1(\\mathbf R))$ norm)")
        plt.legend()
        plt.savefig('comparison-degree-errors.pdf', bbox_inches='tight')
        plt.show()

    except IOError:
        pass

    degrees = list(range(80, degree))
    errors = np.zeros(len(degrees))
    method = args.method if args.method else "ode45"

    for i, d in enumerate(degrees):
        result = solve(method, subdegree=d)
        errors[i] = error_Linf(result)
    np.save("comparison_carrillo_degrees_{}".format(method), degrees)
    np.save("comparison_carrillo_errors_{}".format(method), errors)

# Plots {{{1
if args.test_plots:

    result = solve('ode45')

    times = [0, 1, 2, 5]
    for t in times:
        i = np.argmin(np.abs(np.asarray(time) - t))
        fig, ax = plt.subplots()
        ax.set_xlim(-4, 4)
        quad_comparison.plot(result[i], ax=ax, bounds=False,
                             label="Spectral")
        quad_comparison.plot(rho[i], ax=ax, bounds=False,
                             label="Finite-volume",
                             title="Time: {0:.2f}".format(time[i]))
        ax.set_xlabel('$x$')
        if t == 5:
            plt.legend(loc='lower left')
        plt.savefig('comparison_carrillo_solution-{}.pdf'.format(t),
                    bbox_inches='tight')
        plt.show()

    plt.ion()
    fig, ax = plt.subplots(2)

    def plot(i, t, title=None):

        ax[0].clear()
        quad_comparison.plot(t, ax=ax[0], bounds=False, title=title)
        quad_comparison.plot(rho[i], ax=ax[0], bounds=False, title="Solutions")
        ax[0].set_xlim(-4, 4)

        ax[1].clear()
        evalualion = quad_comparison.eval(t)
        quad_comparison.plot(rho[i] - evalualion, ax=ax[1],
                             bounds=False, title="Error")
        # t.plot(ax=ax[1])

        plt.draw()
        plt.pause(0.01)

    for i, t in enumerate(result):
        Δ = quad_comparison.eval(t) - rho[i]
        Δ = np.sum(np.abs(Δ))*(xs[-1] - xs[0])/(len(Δ) - 1)
        print(Δ)
        if i % 10 == 0:
            plot(i, t)

# Test animations {{{1
if args.test_animate:
    result = solve('ode45')
    writer = animation.writers['ffmpeg'](fps=15, bitrate=1800)
    fig, ax = plt.subplots()
    ax.set_xlabel('$x$')
    l1, = ax.plot([], [], label='Spectral method')
    l2, = ax.plot([], [], label='Finite volume method')
    time_text = ax.text(.08, .05, "Time: 0.00", fontsize=18,
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes)
    legend = ax.legend()

    def init():
        ax.set_xlim(-4, 4)
        ax.set_ylim(-.1, .6)
        return time_text, legend, l1, l2

    fig_tmp, ax_tmp = plt.subplots()
    factor = 10

    def update(i):
        index = i*factor
        ax.set_prop_cycle(None)
        newl1, = quad_comparison.plot(result[index], ax=ax_tmp)
        newl2, = quad_comparison.plot(rho[index], ax=ax_tmp)
        l1.set_data(*newl1.get_data())
        l2.set_data(*newl2.get_data())
        time_text.set_text('Time: {0:.2f}'.format(time[index]))
        return time_text, legend, l1, l2

    anim = animation.FuncAnimation(fig, update, len(result) // factor,
                                     init_func=init, repeat=False, blit=True)
    anim.save('white-noise-d={}.avi'.format(degree), writer=writer)

# Test animation of Hermite coefficients {{{1
if args.test_animate_coefficients:
    result = solve('ode45')
    writer = animation.writers['ffmpeg'](fps=15, bitrate=1800)
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xlabel('d')
    ax.set_ylim(1e-3,1)
    rects = ax.bar(range(degree + 1), abs(result[0].coeffs), bottom=0.001)
    time_text = ax.text(.92, .95, "Time: 0.00", fontsize=18,
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=ax.transAxes)
    fig_tmp, ax_tmp = plt.subplots()
    factor = 10

    def init():
        return rects, time_text
    def update(i):
        index = i*factor
        new_rects = result[index].plot(ax=ax_tmp)
        for r, newr in zip(rects, new_rects):
            r.set_height(abs(newr.get_height()))
        time_text.set_text('Time: {0:.2f}'.format(time[index]))
        return rects, time_text

    anim = animation.FuncAnimation(fig, update, len(result) // factor,
                                   init_func=init, repeat=False, interval=5)
    anim.save('white-noise-coeffs-d={}.avi'.format(degree), writer=writer)
