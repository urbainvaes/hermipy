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
equation = eq.McKean_Vlasov_harmonic_noise
x, y, z, f = equation.x, equation.y, equation.z, equation.f
params = equation.params()
params['γ'], params['θ'] = 0, 0
forward, params = equation.equation(params), equation.params()
β, Vp, θ, m, ε = (params[k] for k in ('β', 'Vp', 'θ', 'm', 'ε'))

# Map the forward operator to a "backward operator"
factor = sym.exp(- sym.Rational(1, 2) * (y*y + z*z))
operator = eq.map_operator(forward, f, factor).expand()

epsilon = params['ε']
L0 = (operator*epsilon**2).expand().subs(epsilon, 0)
L1 = (operator*epsilon - L0/epsilon).expand().subs(epsilon, 0)
L2 = (operator - L1/epsilon - L0/epsilon**2).expand()

# Quadrature used to solve cell problems
degree, nquad, σx, σy = 10, 20, 1, 1
quad_num = quad.Quad.gauss_hermite(nquad, dirs=[1, 2], mean=[0, 0],
                                   cov=[[σx, 0], [0, σy]])

# Discretization in Hermite space
fyz = sym.Function('fyz')(y, z)

# import ipdb; ipdb.set_trace()
op = quad_num.discretize_op(L0.subs(f, fyz), fyz, degree, 2,
                            sparse=False, index_set="triangle")
# }}}
# Expansion of the solution {{{
nterms = 8
zeros = [sym.Integer(0)] * nterms
u, centered, unk = zeros, zeros.copy(), zeros.copy()


def iL0(term):
    t_rhs = quad_num.transform(term, degree)
    solution = la.solve(op.matrix[1:, 1:], t_rhs.coeffs[1:])
    solution = np.array([0, *solution])
    sol_series = series.Series(solution, t_rhs.position, significant=14)
    symbolic = sol_series.to_function().as_format('xyz')
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

    split = func.Function(rhs.doit().expand(), dim=3, allow_sym=True).split()

    for term in split:
        x_part = term[-1] * term[0].as_format('xyz')
        yz_part = term[1].as_format('xyz') * term[2].as_format('xyz')
        t_rhs = quad_num.transform(yz_part, degree)
        centered[i] += round(t_rhs.coeffs[0], 10) * x_part
        u[i] += x_part * iL0(yz_part)
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
for i in range(3, 6):
    print("Equation for u_{}: ".format(i-2))
    sym.pprint(centered[i])

    # The solution is 0
    for j in range(nterms):
        u[j] = u[j].subs(unk[i-2], 0).doit().expand()
        centered[j] = centered[j].subs(unk[i-2], 0).doit().expand()

# Centering condition for u₂
print("Equation for u₄: ")
sym.pprint(centered[6])

C1, C2 = sym.symbols('C1 C2', real=True)

fokker_planck = - ((1/β)*sym.exp(-β*Vp)*(unk[4]*sym.exp(β*Vp)).diff(x)).diff(x)
remainder = centered[6] - fokker_planck

beta3 = (remainder*β**3).expand().subs(β, 0)
beta2 = ((remainder - beta3/β**3) * β**2).expand().subs(β, 0)
beta1 = ((remainder - beta3/β**3 - beta2/β**2) * β).expand().subs(β, 0)
beta0 = remainder - beta3/β**3 - beta2/β**2 - beta1/β

assert (unk[0].diff(x)*Vp.diff(x)).diff(x, x, x, x) \
        + (unk[0].diff(x, x, x)*Vp.diff(x)).diff(x, x) == beta2

assert (Vp.diff(x)*Vp.diff(x)*unk[0].diff(x)).diff(x, x, x) \
        - (Vp.diff(x)*Vp.diff(x, x)*unk[0].diff(x)).diff(x, x) - beta1 == 0

int_fokker_planck = - (1/β)*sym.exp(-β*Vp)*(unk[4]*sym.exp(β*Vp)).diff(x)

int_beta1 = (Vp.diff(x)*Vp.diff(x)*unk[0].diff(x)).diff(x, x) \
            - (Vp.diff(x)*Vp.diff(x, x)*unk[0].diff(x)).diff(x)

int_beta2 = (unk[0].diff(x)*Vp.diff(x)).diff(x, x, x) \
            + (unk[0].diff(x, x, x)*Vp.diff(x)).diff(x)

int_beta3 = unk[0].diff(x, x, x, x, x)

integral1 = int_fokker_planck + \
            int_beta1/β + int_beta2/β**2 + int_beta3/β**3 + C1

assert (integral1.diff(x) - centered[6]).expand() == 0

# C1 = 0
integral1 = integral1.subs(C1, 0)

integrand = ((integral1.subs(unk[0], solution_0).doit().expand()
                       .subs(unk[4], solution_0*unk[4]).doit().expand()
              / solution_0).expand())

integral2 = 3*(Vp.diff(x, x))**2/(2*β) - Vp.diff(x, x, x, x)/β**2 \
            - Vp.diff(x, x)*Vp.diff(x)**2 - unk[4]/β \
            + 2/β*Vp.diff(x)*Vp.diff(x, x, x) - Vp.diff(x, x)**2/β \
            + (Vp.diff(x)*Vp.diff(x, x)**2).integrate(x) + C2

assert (integrand - integral2.diff(x)).expand() == 0

solution_4 = sym.solve(integral2, unk[4])[0] * solution_0

assert centered[6].subs(unk[0], solution_0)\
                  .subs(unk[4], solution_4).doit().expand() == 0
# }}}

solution = 0
for i in range(5):
    solution += epsilon**i * u[i]

print("Solution: ")
sym.pprint(func.Function.sanitize(solution).factor())

solution_x = 0
split = func.Function(solution.expand(), dim=3, allow_sym=True).split()

for term in split:
    x_part = term[-1] * term[0].as_format('xyz')
    yz_part = term[1].as_format('xyz') * term[2].as_format('xyz')
    solution_x += x_part * quad_num.integrate(yz_part)

solution_x = func.Function.sanitize(solution_x)\
            .subs(unk[0], solution_0)\
            .subs(unk[4], solution_4)

print("x-projection of the solution: ")
sym.pprint(func.Function.sanitize(solution_x).factor())

# Checks
solution = solution.subs(unk[0], solution_0).subs(unk[4], solution_4)

for i in range(nterms):
    u[i] = u[i].subs(unk[0], solution_0)\
               .subs(unk[4], solution_4)

assert (operator/epsilon**4).expand()\
        .subs(f, solution + ε**5*u[5] + ε**6*u[6])\
        .doit().expand().subs(epsilon, 0) == 0


# Particular
potential = (x**4/4 - x**2/2) + θ/2*(x - m)**2

θn, mn, βn, εn = .5, sym.Rational(1, 5), sym.Rational(1, 10), .5

solution_x_n = solution_x.subs(Vp, potential).doit().expand()
for i in range(nterms):
    u[i] = u[i].subs(Vp, potential).doit().expand()

solution_x_n = solution_x_n.factor().subs(((θ, θn), (m, mn), (β, βn), (ε, εn)))
for i in range(nterms):
    u[i].subs(((θ, θn), (m, mn), (β, βn), (ε, εn)))

import sympy.plotting as sp
sp.plot(solution_x_n)

quadx = quad.Quad.gauss_hermite(nquad, dirs=[0], mean=[0], cov=[[1]])
solution_p = solution.subs(Vp, potential)
