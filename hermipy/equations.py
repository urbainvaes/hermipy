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

import hermipy.quad as quad
import sympy as sym
import numpy as np
import numpy.linalg as la
sym.init_printing()


def map_operator(operator, function, factor):
    return (operator.subs(function, function*factor).doit()/factor).expand()


def solve_gaussian(operator, f, variables):

    if len(variables) is 3:
        x, y, z = variables
        cx, cy, cz = sym.symbols('cx cy cz', real=True)
        cxy, cyz, cxz = sym.symbols('cxy cyz cxz', real=True)
        ansatz = sym.exp(- sym.Rational(1, 2) *
                         (cx*x*x + cz*z*z + cy*y*y
                          + 2*cxy*x*y + 2*cyz*y*z + 2*cxz*x*z))
        res = (operator.subs(f, ansatz)/ansatz).doit().cancel()
        coeffs_dict = res.as_poly(x, y, z).as_dict()
        coeffs_list = list(coeffs_dict[k] for k in coeffs_dict)
        (sol_cx, sol_cy, sol_cz, sol_cxy, sol_cyz, sol_cxz) = \
            sym.solve(coeffs_list, cx, cy, cz, cxy, cyz, cxz)[0]
        solution = ansatz.subs([(cx, sol_cx), (cy, sol_cy),
                                (cz, sol_cz), (cxy, sol_cxy),
                                (cyz, sol_cyz), (cxz, sol_cxz)])
        determinant = sol_cx*(sol_cy*sol_cz-sol_cyz**2) \
            + sol_cxy*(sol_cxz*sol_cyz-sol_cxy*sol_cz) \
            + sol_cxz*(sol_cxy*sol_cyz-sol_cy*sol_cxz)
        factor = sym.sqrt(determinant / (2*sym.pi)**3)
        solution = factor * solution

    if len(variables) is 2:
        x, y = variables
        cx, cxy, cy = sym.symbols('sol_cx sol_cxy sol_cy', real=True)
        ansatz = sym.exp(-sym.Rational(1, 2)*(cx*x*x + 2*cxy*x*y + cy*y*y))
        res = (operator.subs(f, ansatz)/ansatz).doit().cancel()
        coeffs_dict = res.as_poly(x, y).as_dict()
        coeffs_list = list(coeffs_dict[k] for k in coeffs_dict)
        (sol_cx, sol_cxy, sol_cy) = sym.solve(coeffs_list, cx, cxy, cy)[0]
        solution = ansatz.subs([(cx, sol_cx), (cxy, sol_cxy), (cy, sol_cy)])
        determinant = sol_cx * sol_cy - sol_cxy * sol_cxy
        factor = sym.sqrt(determinant) / (2*sym.pi)
        solution = factor * solution

    if len(variables) is 1:
        x = variables[0]
        s2 = sym.symbols('s2', real=True, positive=True)
        ansatz = 1/sym.sqrt(2*sym.pi*s2)*sym.exp(-sym.Rational(1, 2)*(x*x/s2))
        res = (operator.subs(f, ansatz)/ansatz).doit().cancel()
        coeffs_dict = res.as_poly(x).as_dict()
        coeffs_list = list(coeffs_dict[k] for k in coeffs_dict)
        sol_s2 = sym.solve(coeffs_list, s2)[s2]
        solution = ansatz.subs(s2, sol_s2)

    assert operator.subs(f, solution).doit().expand().cancel() == 0
    return solution


class Fokker_Planck_1d:

    # Space variables
    variables = sym.symbols('x', real=True)

    # Sharthand notations
    x = variables

    # Unknown
    f = sym.Function('f')(x)

    @classmethod
    def params(cls):
        return {'β': sym.symbols('β', real=True, positive=True),
                'Vp': sym.Function('Vp')(cls.x)}

    @classmethod
    def equation(cls, params):

        β, Vp = (params[x] for x in ['β', 'Vp'])

        # Shorthand notations
        d, x, f = sym.diff, cls.x, cls.f

        # Fokker planck operator
        return d(d(Vp, x)*f + 1/β * d(f, x), x)


class McKean_Vlasov_white:

    # Space variables
    x = sym.symbols('x', real=True)

    # Unknown
    f = sym.Function('f')(x)

    @classmethod
    def params(cls):

        params, functions = {}, {}

        # Real positive parameters
        options = {'real': True, 'positive': True}
        for param in ['β', 'γ', 'ε']:
            params[param] = sym.symbols(param, **options)

        # Real parameters
        options = {'real': True}
        for param in ['θ', 'm']:
            params[param] = sym.symbols(param, **options)

        # Potential functions
        functions['Vp'] = sym.Function('Vp')(cls.x)

        return {**params, **functions}

    @classmethod
    def fluxes(cls, params):

        # Shorthand notations
        d, x, f = sym.diff, cls.x, cls.f

        # Real parameters
        β, ε, θ, m = (params[x] for x in ['β', 'ε', 'θ', 'm'])

        # Functional parameter
        Vp = params['Vp']

        # Fokker planck operator
        return - (d(Vp, x)*f + θ*(x-m)*f + 1/β * d(f, x))

    @classmethod
    def equation(cls, params):

        # Compute fluxes
        flux = cls.fluxes(params)

        # Fokker planck operator
        return - flux.diff(cls.x)


class McKean_Vlasov:

    # Space variables
    variables = sym.symbols('x y', real=True)

    # Sharthand notations
    x, y = variables

    # Unknown
    f = sym.Function('f')(variables[0], variables[1])

    @classmethod
    def params(cls):

        params, functions = {}, {}

        # Real positive parameters
        options = {'real': True, 'positive': True}
        for param in ['β', 'γ', 'ε']:
            params[param] = sym.symbols(param, **options)

        # Real parameters
        options = {'real': True}
        for param in ['θ', 'm']:
            params[param] = sym.symbols(param, **options)

        # Potential functions
        functions['Vp'] = sym.Function('Vp')(cls.x)
        functions['Vy'] = sym.Function('Vy')(cls.y)

        return {**params, **functions}

    @classmethod
    def fluxes(cls, params):

        # Shorthand notations
        d, x, y, f = sym.diff, cls.x, cls.y, cls.f

        # Real parameters
        β, γ, ε, θ, m = (params[x] for x in ['β', 'γ', 'ε', 'θ', 'm'])

        # Functional parameter
        Vp, Vy = params['Vp'], params['Vy']

        # Effective diffusion
        functions = Vy.atoms(sym.function.AppliedUndef)

        if functions == set():
            degree, n_points, q = 100, 201, quad.Quad.gauss_hermite
            fy = sym.Function('f')(y)
            gen = Vy.diff(y)*fy.diff(y) - fy.diff(y, y)
            σy, density = .2, sym.exp(-Vy)
            gaussian = sym.exp(-y*y/(2*σy))/sym.sqrt(2*sym.pi*σy)
            integral = q(n_points, dirs=[1]).integrate(density, flat=True)
            factor = sym.sqrt(gaussian / (density / integral))
            qy = q(n_points, dirs=[1], factor=factor, cov=[[σy]])
            L0 = qy.discretize_op(gen, degree, sparse=False,
                                  index_set="triangle")
            rhs = qy.transform('y', degree)
            solution = la.solve(L0.matrix[1:, 1:], rhs.coeffs[1:])
            coeff_noise = np.dot(solution, rhs.coeffs[1:])
            coeff_noise = sym.Rational(coeff_noise).limit_denominator(1e8)
            print("Effective noise: " + str(float(coeff_noise)))
            coeff_noise = sym.sqrt(1/β/coeff_noise)

        else:
            # Assume Vy is y*y/2
            coeff_noise = sym.sqrt(1/β)

        # Fokker planck operator
        flux_x = - (d(Vp, x)*f + θ*(x-m)*f - (1-γ)*coeff_noise*y*f/ε
                    + γ**2/β * d(f, x))
        flux_y = - (1/ε**2) * (f*Vy.diff(y) + d(f, y))
        return [flux_x, flux_y]

    @classmethod
    def equation(cls, params):

        # Shorthand notations
        x, y = cls.x, cls.y

        # Compute fluxes
        fx, fy = cls.fluxes(params)

        # Fokker planck operator
        return - fx.diff(x) - fy.diff(y)


class McKean_Vlasov_harmonic_noise:

    # Space variables
    variables = sym.symbols('x y z', real=True)

    # Sharthand notations
    x, y, z = variables

    # Unknown
    f = sym.Function('f')(x, y, z)

    @classmethod
    def params(cls):

        params, functions = {}, {}

        # Real positive parameters
        options = {'real': True, 'positive': True}
        for param in ['β', 'γ', 'ε']:
            params[param] = sym.symbols(param, **options)

        # Real parameters
        options = {'real': True}
        for param in ['θ', 'm']:
            params[param] = sym.symbols(param, **options)

        # Potential function
        functions['Vp'] = sym.Function('Vp')(cls.x)

        return {**params, **functions}

    @classmethod
    def fluxes(cls, params):

        # Real parameters
        β, γ, ε, θ, m = (params[x] for x in ['β', 'γ', 'ε', 'θ', 'm'])

        # Functional parameter
        Vp = params['Vp']

        # Shorthand notations
        d, x, y, z, f = sym.diff, cls.x, cls.y, cls.z, cls.f

        # Friction coefficient
        λ = 1

        # Fokker planck operator
        flux_x = - (d(Vp, x)*f + θ*(x-m)*f - (1-γ)*sym.sqrt(1/β)*z*f/ε
                    + γ**2/β * d(f, x))
        flux_y = - (1/ε**2) * (f*z + λ*(f*y + d(f, y)))
        flux_z = - (1/ε**2) * (-y*f)
        return [flux_x, flux_y, flux_z]

    @classmethod
    def equation(cls, params):

        # Shorthand notations
        x, y, z = cls.x, cls.y, cls.z

        # Compute fluxes
        fx, fy, fz = cls.fluxes(params)

        # Fokker planck operator
        return - fx.diff(x) - fy.diff(y) - fz.diff(z)


class Langevin:

    # Space variables
    variables = sym.symbols('x y', real=True)

    # Sharthand notations
    x, y = variables

    # Unknown
    f = sym.Function('f')(x, y)

    @classmethod
    def forward(cls, λ):

        # Shorthand notations
        d, x, y, f = sym.diff, cls.x, cls.y, cls.f

        # Fokker planck operator
        operator = (-x*d(f, y) + y*d(f, x)) + λ*d(f*x + d(f, x), x)

        return operator

    @classmethod
    def backward(cls, λ):

        # Shorthand notations
        d, x, y, f = sym.diff, cls.x, cls.y, cls.f

        # Backward Kolmogorov operator
        operator = (x*d(f, y) - y*d(f, x)) + λ*(- x*d(f, x) + d(f, x, x))

        return operator

# vim: foldmethod=indent
