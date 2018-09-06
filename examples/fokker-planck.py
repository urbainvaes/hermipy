#1 Copyright (C) 2018 Urbain Vaes

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
parser.add_argument('-dir', '--directory', type=str)
parser.add_argument('-i', '--interactive', action='store_false')
parser.add_argument('-tc', '--convergence', action='store_true')
parser.add_argument('-b', '--beta', type=float)
parser.add_argument('-l', '--lambda_f', type=float)
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

# Dimension of the problem
dim = 1

# Equation parameters
r, s = sym.Rational, sym.symbols
equation = eq.Fokker_Planck_1d
x, f = equation.x, equation.f


def set_param(arg, default, symbol):
    if arg and r(arg) < 0:
        return s(symbol)
    if arg:
        return r(arg)
    return default


# Fokker Planck equation
β = set_param(args.beta, r(1), 'β')
λ = set_param(args.lambda_f, r(1, 2), 'λ')
Vp = sym.Piecewise((1 + 2*sym.log(-x), x < -1),
                   (x**2, x < 1), (1 + 2*sym.log(x), True))
params = {'β': β, 'Vp': Vp}
forward = equation.equation(params)

sx, mx = r(1, 10), r(0)
degree = args.degree if args.degree else 100
kwargs0 = {'degree': degree}
n_points_num = 2*degree + 1

# Potential for approximation
Vq = sym.Rational(1/2)*x*x/sx
# factor = sym.exp(-(λ*Vq + β*(1-λ)*Vp))
factor = sym.exp(-Vq/2)

# Calculation of the solution
new_q = hm.Quad.gauss_hermite
mean, cov = [mx], [[sx]]
quad = new_q(n_points_num, mean=mean, cov=cov, factor=factor)

# Integral and moment operator
w = quad.factor * quad.factor / quad.position.weight()
In = quad.transform(w, degree=degree)
m1 = quad.transform(w * x, degree=degree)
m2 = quad.transform(w * x * x, degree=degree)

# Initial condition
u = sym.exp(-x*x/2)
u = u / quad.integrate(u, flat=True)

# Degrees
degrees = range(10, degree)
varf = quad.discretize_op(forward, degree=degree)
s_exact = sym.exp(-β*Vp)
t_exact = quad.transform(s_exact, degree=degree)
integral = float(In*t_exact)
s_exact, t_exact = s_exact / integral, t_exact / integral
eig_ground, errors, mins = [], [], []
if args.interactive:
    quad.plot(s_exact)

# Initial condition
t = quad.transform(u, degree=degree)

if args.interactive:
    plt.ion()
    fig, ax = plt.subplots(2)


def plot(t, title=None):

    if not args.interactive:
        return

    ax[0].clear()
    quad.plot(t_exact, ax=ax[0], bounds=False)
    quad.plot(t, ax=ax[0], bounds=True, title=title)

    ax[1].clear()
    t.plot(ax=ax[1])

    plt.draw()
    plt.pause(0.01)


eye = quad.varf('1', degree=degree)
t = t_exact
for d in degrees:
    sub_varf = varf.subdegree(d)
    sub_eye = eye.subdegree(d)
    Id = In.subdegree(d)
    t = t.subdegree(d)

    dt, Ns, Δ = 2e-4, int(1e4), np.inf
    for i in range(Ns):

        if args.interactive:
            plot(t)

        print("d: " + str(d) + ", i: " + str(i) +
              ", dt: " + str(dt) + ", Δ: " + str(Δ))

        # Backward Euler
        total_op = sub_eye - dt*sub_varf
        new_t = total_op.solve(t)

        # Normalization
        new_t = new_t / float(Id*new_t)

        # Error
        Δ = quad.norm(new_t - t, n=1, flat=True)

        # Time adaptation
        threshold, dt_max = .01, 64
        if Δ*dt > 2*threshold:
            dt = dt / 2.
            continue
        elif Δ*dt < threshold and dt < dt_max:
            dt = dt * 2.
        t = new_t

        if Δ < 1e-10:

            error = quad.norm(t.subdegree(degree) - t_exact, n=1, flat=True)
            min = np.min(quad.eval(t))
            eig = quad.norm(t-sub_varf(t), n=1, flat=True) \
                / quad.norm(t, n=1, flat=True)
            mins.append(abs(min))
            errors.append(error)
            print(min, error)

            if args.interactive:
                plot(t)

                d = d - 1
                break


# for d in degrees:
#     print("-- Solving for degree = " + str(d))
#     sub_varf = varf.subdegree(d)
#     [l], [e] = sub_varf.eigs(k=1, which='SM')

#     # Normalize and ensure sign consistency
#     e = e / float(In.subdegree(d)*e)
#     error = quad.norm(e.subdegree(degree) - e_exact, n=1, flat=True)
#     plot(e, title=str(error))

#     mins.append(np.min(quad.eval(e)))
#     eig_ground.append(l)
#     errors.append(error)
#     print(error, l, np.min(quad.eval(e)))
