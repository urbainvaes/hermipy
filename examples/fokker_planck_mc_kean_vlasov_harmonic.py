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

import ipdb
import sympy as sym
import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as las
import matplotlib.pyplot as plt

import hermipy.quad as hm
import hermipy.equations as eq
import hermipy.settings as rc
import hermipy.core as core

rc.settings['tensorize'] = True
rc.settings['trails'] = True

dim = 3

# Equation parameters
r = sym.Rational
equation = eq.McKean_Vlasov_harmonic_noise
x, y, z, f = equation.x, equation.y, equation.z, equation.f
params = {'β': r(10), 'ε': r(.5), 'γ': 0, 'θ': 0, 'm': 0}
Vp, degree, β = x**4/4 - x**2/2, 10, params['β']

# Numerical parameters
s2x, s2y, s2z, degree = r(1, 10), 1, 1, 15
n_points_num = 2*degree + 1

# Calculation of the solution
new_q = hm.Quad.gauss_hermite
cov = [[s2x, 0, 0], [0, s2y, 0], [0, 0, s2z]]
quad = new_q(n_points_num, dim=3, mean=[0]*3, cov=cov)

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

# Mapped operator
backward = eq.map_operator(forward, f, factor)

# Discretize factor
factor_x = quad.project(0).discretize(factor_x)
factor = quad.discretize(factor)

degrees = list(range(5, degree))

# Discretization of the operator
mat = quad.discretize_op(backward, degrees[-1], sparse=True)

solutions = []

v0, eig_vec = None, None
for d in degrees:
    print(d)
    npolys = core.iterator_size(dim, degree)
    if d is not degrees[0]:
        v0 = np.zeros(npolys)
        for i in range(len(eig_vec)):
            v0[i] = eig_vec[i]
    sub_mat = (mat[0:npolys, 0:npolys]).copy(order='C')
    # pdb.set_trace()
    eig_vals, eig_vecs = las.eigs(sub_mat, k=1, v0=v0, which='LR')
    eig_vec = np.real(eig_vecs.T[0])
    ground_state = eig_vec * np.sign(eig_vec[0])
    ground_state_eval = quad.eval(quad.series(ground_state))*factor
    norm = quad.norm(ground_state_eval, n=1, flat=True)
    ground_state_eval = ground_state_eval / norm
    solutions.append(ground_state_eval)

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
