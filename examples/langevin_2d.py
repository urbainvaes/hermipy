# Import necessary modules
import math
import time
import sympy as sym
import numpy as np
import hermipy as hm
import scipy.integrate as integrate

# Set Hermipy settings
hm.settings['tensorize'] = True
hm.settings['sparse'] = True
hm.settings['cache'] = False

# Variables, corresponding to (q, p, z1, z2), and unknown function
x, y, z, w = sym.symbols('x y z w')
f = sym.Function('f')(x, y, z, w)

# Vector containing auxiliary variables
z_vec = sym.Matrix([z, w])

# Physical parameters
b = 1
g = 1
d = 1e-10
# d = 0

# Periodic potential
delta = sym.symbols('δ')
V0 = (1/2)*(1 - sym.cos(x)) + (1/2)*(1 - sym.cos(y))
V1 = (1/2)*sym.cos(x)*sym.cos(y)
V = V0 + delta*V1

# V = (1/2)*(1 - sym.cos(y))
# V = (1/2)*(1 - sym.cos(x))

# Normalization
lambda_V = sym.lambdify((x, y), sym.exp(-b*(V0 + d*V1)))
zx = integrate.dblquad(lambda_V, -math.pi, math.pi,
                       lambda x: -math.pi, lambda x: math.pi)[0]
print("Normalization: {}".format(zx))

# Backward Kolmogorov operator
backward = \
    z*f.diff(x) + w*f.diff(y) \
    - V.diff(x)*f.diff(z) - V.diff(y)*f.diff(w) \
    + g*(- z*f.diff(z) + (1/b)*f.diff(z, z)) \
    + g*(- w*f.diff(w) + (1/b)*f.diff(w, w))


L0 = backward.subs(delta, 0)
L1 = ((backward - L0)/delta).cancel()

# Set parameters for quadrature
kwargs = {'degree': 2, 'index_set': 'cube'}
npoints = kwargs['degree']*2 + 1

# # Define quadratures
# qxy = hm.Quad.fourier(npoints, dirs=[0, 1],
#                       factor=math.sqrt(zx) * sym.exp(b/2*(V0 + d*V1)))
# qz = hm.Quad.gauss_hermite(npoints, dirs=[2], cov=[[1]], factor=1)
# qw = hm.Quad.gauss_hermite(npoints, dirs=[3], cov=[[1]], factor=1)

# # Calculate the tensor product of the quadratures
# # quad = qx*qy*qz*qw
# quad = qxy*qz*qw

# # Calculate the matrix discretizing 'backward'
# M0 = quad.discretize_op(L0, **kwargs)
# M1 = quad.discretize_op(L1, **kwargs)
# discretized_op = M0 + d*M1

# # Calculate the right-hand side
# rhs1 = quad.transform('z', **kwargs)
# rhs2 = quad.transform('w', **kwargs)

# # Calculate the projection of the function in the kernel of 'backward'
# remove = quad.transform(1, **kwargs)

# # Calculate solution using the saddle-point formulation
# solution1 = (-discretized_op).solve(rhs1, remove_vec=remove)
# solution2 = (-discretized_op).solve(rhs2, remove_vec=remove)

# # Calculate the diffusion coefficient
# diffusion_11 = float(solution1*rhs1)
# diffusion_12 = float(solution1*rhs2)
# diffusion_21 = float(solution2*rhs1)
# diffusion_22 = float(solution2*rhs2)
# diffusion_matrix = np.matrix([[diffusion_11, diffusion_12],
#                               [diffusion_21, diffusion_22]])
# print(diffusion_matrix)
# time.sleep(1)

f2 = sym.Function('f')(x, y)
fx = sym.Function('f')(x)
kwargs = {'degree': 2, 'index_set': 'cube'}
operator = sym.cos(x)*f2
hm.settings['tensorize'] = True
qx = hm.Quad.fourier(npoints, dirs=[0], factor=1)
qy = hm.Quad.fourier(npoints, dirs=[1], factor=1)
# qxy_1 = hm.Quad.fourier(npoints, dirs=[0, 1],
#                         factor=math.sqrt(zx) * sym.exp(b/2*(V0 + 1e-10*V1)))
# qxy_2 = hm.Quad.fourier(npoints, dirs=[0, 1],
#                         factor=math.sqrt(zx) * sym.exp(b/2*(V0)))
# hm.settings['tensorize'] = True
# o1 = (qx*qy).discretize_op(operator, **kwargs)

# hm.settings['tensorize'] = False
# o2 = (qx*qy).discretize_op(operator, **kwargs)

# print(o2.list

# print(qx.discretize_op(sym.sin(x)*fx, **kwargs).matrix.todense())
# print(o1.matrix.todense())
# print(o2.matrix.todense())
# print(np.max(np.abs(o1.matrix - o2.matrix)))

print(((qx*qy).transform('sin(x)', tensorize=False, **kwargs)).coeffs)
print(((qx*qy).varf('sin(x)', tensorize=False, **kwargs)).matrix.todense())
# print(((qx*qy).varf('sin(x)', tensorize=True, **kwargs)).matrix.todense())
