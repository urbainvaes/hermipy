# Import necessary modules
import math
import argparse
import sympy as sym
import numpy as np
import hermipy as hm
import scipy.integrate as integrate

parser = argparse.ArgumentParser()
parser.add_argument('-td', '--test_delta', action='store_true')
parser.add_argument('-tg', '--test_gamma', action='store_true')
parser.add_argument('-g', '--gamma', type=str)
parser.add_argument('-d', '--delta', type=str)
args = parser.parse_args()

# Set Hermipy settings
hm.settings['tensorize'] = True
hm.settings['sparse'] = True
hm.settings['cache'] = True

# Variables, corresponding to (q, p, z1, z2), and unknown function
x, y, z, w = sym.symbols('x y z w')
f = sym.Function('f')(x, y, z, w)
fxy = sym.Function('f')(x, y)

# Vector containing auxiliary variables
z_vec = sym.Matrix([z, w])

# Physical parameters
b = 1
gamma = sym.Symbol('γ')
delta = sym.Symbol('δ')

# Periodic potential
V0 = (1/2)*(1 - sym.cos(x)) + (1/2)*(1 - sym.cos(y))
V1 = (1/2)*sym.cos(x)*sym.cos(y)
V = V0 + delta*V1

# Backward Kolmogorov operator
backward = \
    z*f.diff(x) + w*f.diff(y) \
    - V.diff(x)*f.diff(z) - V.diff(y)*f.diff(w) \
    + gamma*(- z*f.diff(z) + (1/b)*f.diff(z, z)) \
    + gamma*(- w*f.diff(w) + (1/b)*f.diff(w, w))

# Set parameters for quadrature
kwargs = {'degree': 8, 'index_set': 'cube'}
npoints = kwargs['degree']*2 + 1


def overdamped(d):
    lambda_V = sym.lambdify((x, y), sym.exp(-b*(V0 + d*V1)))
    zx = integrate.dblquad(lambda_V, -math.pi, math.pi,
                           lambda x: -math.pi, lambda x: math.pi)[0]
    qxy = hm.Quad.fourier(npoints, dirs=[0, 1],
                          factor=math.sqrt(zx) * sym.exp(b/2*(V0 + d*V1)))
    backward_overdamped = \
        - V.diff(x)*fxy.diff(x) - V.diff(y)*fxy.diff(y) \
        + (1/b)*fxy.diff(x, x) + (1/b)*fxy.diff(y, y)

    Mo = qxy.discretize_op(backward_overdamped.subs(delta, d), **kwargs)
    rhs1 = qxy.transform(-V.subs(delta, d).diff(x), **kwargs)
    rhs2 = qxy.transform(-V.subs(delta, d).diff(y), **kwargs)
    remove = qxy.transform(1, **kwargs)
    solution1 = (-Mo).solve(rhs1, remove_vec=remove)
    solution2 = (-Mo).solve(rhs2, remove_vec=remove)
    Iop = qxy.transform(1, **kwargs)
    dx = qxy.discretize_op(fxy.diff(x), **kwargs)
    dy = qxy.discretize_op(fxy.diff(y), **kwargs)
    diffusion1_o = ((1/b) + float(solution1*rhs1) +
                    (2/b)*float(Iop*dx(solution1)))
    diffusion2_o = ((1/b) + float(solution2*rhs2) +
                    (2/b)*float(Iop*dy(solution2)))
    return diffusion1_o, diffusion2_o


def diffusion(d, g):
    operator = backward.subs(((delta, d), (gamma, g)))

    # Normalization
    lambda_V = sym.lambdify((x, y), sym.exp(-b*(V0 + d*V1)))
    zx = integrate.dblquad(lambda_V, -math.pi, math.pi,
                           lambda x: -math.pi, lambda x: math.pi)[0]

    # Define quadratures
    qz = hm.Quad.gauss_hermite(npoints, dirs=[2], cov=[[1]], factor=1)
    qw = hm.Quad.gauss_hermite(npoints, dirs=[3], cov=[[1]], factor=1)
    qxy = hm.Quad.fourier(npoints, dirs=[0, 1],
                          factor=math.sqrt(zx) * sym.exp(b/2*(V0 + d*V1)))

    # Calculate the tensor product of the quadratures
    quad = qxy*qz*qw

    # Calculate the matrix discretizing 'backward'
    discretized_op = quad.discretize_op(operator, **kwargs)

    # Calculate the right-hand side
    rhs1 = quad.transform('z', **kwargs)
    rhs2 = quad.transform('w', **kwargs)

    # Calculate the projection of the function in the kernel of 'backward'
    remove = quad.transform(1, **kwargs)

    # Calculate solution using the saddle-point formulation
    solution1 = (-discretized_op).solve(rhs1, remove_vec=remove)
    solution2 = (-discretized_op).solve(rhs2, remove_vec=remove)

    # Calculate the diffusion coefficient
    diffusion_11 = float(solution1*rhs1)
    diffusion_12 = float(solution1*rhs2)
    diffusion_21 = float(solution2*rhs1)
    diffusion_22 = float(solution2*rhs2)
    diffusion_matrix = np.matrix([[diffusion_11, diffusion_12],
                                  [diffusion_21, diffusion_22]])

    return diffusion_matrix


if args.test_delta:
    asym = 0.815
    g = float(args.gamma)
    ds = np.arange(.005, 1, .05)
    Ds = 0*ds
    for i, d in enumerate(ds):
        Ds[i] = diffusion(d, g)[0, 0]
        print("δ = {}, D = {}".format(ds[i], Ds[i]))

if args.test_gamma:
    d = float(args.delta)
    gs = np.logspace(0, 3, 20)
    Ds, Do = 0*gs, overdamped(d)
    for i, g in enumerate(gs):
        Ds[i] = diffusion(d, g)[0, 0]
        print("γ = {}, D = {}".format(gs[i], Ds[i]))
