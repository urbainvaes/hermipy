#!/usr/bin/env python
#
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

import hermipy.equations as eq
import hermipy.quad as quad
import hermipy.series as series
import hermipy.function as func

import numpy as np
import numpy.linalg as la
import sympy as sym

# Declare variables, operators and parameters {{{
sym.init_printing()
equation = eq.McKean_Vlasov
x, y, f = equation.x, equation.y, equation.f
params = equation.params()
params['γ'], params['θ'] = 0, 0
params.update({'γ': 0, 'θ': 0, 'Vy': y*y/2})
forward = equation.equation(params)
β, Vp, θ, m = (params[k] for k in ('β', 'Vp', 'θ', 'm'))

# Map the forward operator to a "backward operator"
factor = sym.exp(- sym.Rational(1, 2) * (y*y))
operator = eq.map_operator(forward, f, factor).expand()

epsilon = params['ε']
L0 = (operator*epsilon**2).expand().subs(epsilon, 0)
L1 = (operator*epsilon - L0/epsilon).expand().subs(epsilon, 0)
L2 = (operator - L1/epsilon - L0/epsilon**2).expand()

# Quadrature used to solve cell problems
degree, nquad, σy = 10, 20, 1
quad_num = quad.Quad.gauss_hermite(nquad, dirs=[1], mean=[0], cov=[[σy]])

# Discretization in Hermite space
fy = sym.Function('fy')(y)

# import ipdb; ipdb.set_trace()
op = quad_num.discretize_op(L0.subs(f, fy), degree, sparse=False)

# }}}
# Expansion of the solution {{{
nterms = 7
zeros = [sym.Integer(0)] * nterms
u, centered, unk = zeros, zeros.copy(), zeros.copy()


def iL0(term):
    t_rhs = quad_num.transform(term, degree)
    solution = la.solve(op.matrix[1:, 1:], t_rhs.coeffs[1:])
    solution = np.array([0, *solution])
    sol_series = series.Series(solution, t_rhs.position, significant=14)
    symbolic = sol_series.to_function().as_xyz()
    symbolic = func.Function.sanitize(symbolic, max_denom=1e8)

    print("--> Solving cell problem with rhs: " + str(term))
    print("----> Solution: " + str(symbolic))
    diff = (L0.subs(f, symbolic).doit() - term)
    if len(diff.expand().free_symbols) > 0:
        print("----> Error: does not match!")
        sym.pprint(diff)
        import ipdb
        ipdb.set_trace()
    return symbolic


for i in range(nterms):
    print("Solving for power " + str(i))

    unk[i] = sym.Function('u{}'.format(i))(x)
    u[i] = unk[i]

    rhs = sym.Integer(0)
    if i > 0:
        rhs += - L1.subs(f, u[i - 1])
    if i > 1:
        rhs += - L2.subs(f, u[i - 2])

    split = func.Function(rhs.doit().expand(), dim=2).split()

    for term in split:
        x_part = term[-1] * term[0].as_xyz()
        y_part = term[1].as_xyz()
        t_rhs = quad_num.transform(y_part, degree)
        centered[i] += round(t_rhs.coeffs[0], 10) * x_part
        u[i] += x_part * iL0(y_part)
    u[i] = func.Function.sanitize(u[i])
    centered[i] = func.Function.sanitize(centered[i])

# }}}
# Some manual calculations {{{
# Operator
fx = sym.Function('fx')(x)
LFP = ((1/β)*sym.exp(-β*Vp)*(fx*sym.exp(β*Vp)).diff(x)).diff(x)

# Centering condition for u₀
print("Equation for u₀: ")
sym.pprint(centered[2])
Z = sym.symbols('Z', real=True)
solution_0 = sym.exp(-β*Vp)/Z
assert centered[2].subs(unk[0], solution_0).doit().cancel() == 0

# Centering condition for u₁
print("Equation for u₁: ")
sym.pprint(centered[3])

# The solution is 0
for i in range(nterms):
    u[i] = u[i].subs(unk[1], 0).doit().expand()
    centered[i] = centered[i].subs(unk[1], 0).doit().expand()

# Centering condition for u₂
print("Equation for u₂: ")
sym.pprint(centered[4])

C1, C2 = sym.symbols('C1 C2')

integral1 = (
        - (1/β)*sym.exp(-β*Vp)*(unk[2]*sym.exp(β*Vp)).diff(x)
        - (1/β)*(Vp.diff(x)*unk[0].diff(x)).diff(x)
        - (1/β**2) * unk[0].diff(x, x, x)
        + C1)

assert (integral1.diff(x) - centered[4]).doit().expand() == 0

integral2 = ((integral1.subs(unk[0], solution_0).doit().expand()
                       .subs(unk[2], solution_0*unk[2]).doit().expand()
             / solution_0).expand().integrate(x).doit() + C2)

solution_2 = sym.solve(integral2, unk[2])[0] * solution_0

assert centered[4].subs(unk[0], solution_0)\
                  .subs(unk[2], solution_2).doit().expand() == 0

# C1 must be 0
solution_2 = solution_2.subs(C1, 0)

for i in range(nterms):
    u[i] = u[i].subs(unk[0], solution_0)\
               .subs(unk[2], solution_2).doit().expand().cancel()

    centered[i] = centered[i].subs(unk[0], solution_0)\
                             .subs(unk[2], solution_2).doit().expand().cancel()

solution = (u[0] + epsilon**2 * u[2])
print("Solution: ")
sym.pprint(func.Function.sanitize(solution).factor())

solution_x = 0
split = func.Function(solution.expand(), dim=2).split()

for term in split:
    x_part = term[-1] * term[0].as_xyz()
    y_part = term[1].as_xyz()
    solution_x += x_part * quad_num.integrate(y_part)

print("x-projection of the solution: ")
sym.pprint(func.Function.sanitize(solution_x).factor())
import ipdb; ipdb.set_trace()
# }}}
# Plot x - y {{{
n_points = 100
potential = (x**4/4 - x**2/2) + θ/2*(x**2 - 2*x*m)
quad_gauss = quad.Quad.gauss_hermite(n_points, dirs=[0, 1],
                                     cov=[[1, 0], [0, 1]])
quadx = quad.Quad.gauss_hermite(100, dirs=[0], cov=[[1]])

r, ε = sym.Rational, epsilon
θn, mn, βn, εn = r(0), r(0, 5), r(1), r(1, 2)
potential_n = potential.subs(((θ, θn), (m, mn), (β, βn), (ε, εn)))
Z_n = quadx.integrate(sym.exp(-βn*potential_n), flat=True)
solution_0_n = sym.exp(-βn*potential_n)/Z_n
solution_2_n = solution_2.subs(Vp, potential).doit()
solution_2_n = solution_2_n.factor().subs(((θ, θn), (m, mn), (β, βn),
                                           (ε, εn), (Z, Z_n)))
C2_n = quadx.integrate(sym.solve(solution_2_n, C2)[0] * solution_0_n, flat=True)
solution_4_n = solution_2_n.subs(C2, C2_n)
assert abs(quadx.integrate(solution_0_n, flat=True) - 1) < 1e-8
assert abs(quadx.integrate(solution_4_n, flat=True) - 0) < 1e-8

rho_z = 1/sym.sqrt(2*sym.pi) * sym.exp(-y*y/2)

un = [0]*3
for i, p in enumerate(un):
    un[i] = (u[i]*rho_z).subs(unk[0], solution_0)\
                        .subs(unk[2], solution_2).doit()
    un[i] = un[i].subs(Vp, potential).doit()
    un[i] = un[i].factor().subs(((θ, θn), (m, mn), (β, βn),
                                 (ε, εn), (Z, Z_n), (C2, C2_n)))
    assert abs(quad_gauss.integrate(un[i], flat=True) - i == 0) < 1e-8


if False:
    import matplotlib
    import matplotlib.pyplot as plt
    matplotlib.rc('font', size=14)
    matplotlib.rc('font', family='serif')
    matplotlib.rc('text', usetex=True)

    fig, axes = plt.subplots(1, 2)
    quad_visu = quad.Quad.newton_cotes([n_points, n_points],
                                       [2, 2], dirs=[0, 1])
    proj_total = 0

    for i in (0, 1):
        cont = quad_visu.plot(un[i+1], ax=axes[i])
        plt.colorbar(cont, ax=axes[i], pad=.01)
        for c in cont.collections:
            c.set_edgecolor("face")
    plt.show()
