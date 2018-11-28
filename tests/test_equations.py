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

import unittest
import sympy as sym
import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as las
import matplotlib.pyplot as plt

import hermipy as hm
import hermipy.core as core
import hermipy.equations as eq
import hermipy.stats as stats
import hermipy.cache as cache


class TestConvergenceFokkerPlanck1d(unittest.TestCase):

    settings = {
        'cache': False,
        'cachedir': '/tmp/test_hermite',
        'tensorize': False,
        'sparse': False,
        'trails': False,
        }

    def setUp(self):

        # Equation parameters
        equation = eq.Fokker_Planck_1d
        self.x, self.f = equation.x, equation.f

        # Inverse temperature
        self.β = sym.Rational(3)

        # Common parameters of numerical approximation
        self.degree = 100
        self.n_points_num = 2*self.degree + 1

        # Reset default settings
        hm.settings.update(self.settings)

    def sym_calc(self, Vp, m, s2):

        # Equation parameters
        equation = eq.Fokker_Planck_1d
        x, f = equation.x, equation.f

        # Calculation of the solution
        new_q = hm.Quad.gauss_hermite

        quad = new_q(self.n_points_num, dim=1, mean=[m], cov=[[s2]])

        # Potential for approximation
        Vq = sym.Rational(1, 2)*(x-m)*(x-m)/(self.β*s2)

        # Fokker Planck operator
        forward = equation.equation({'β': self.β, 'Vp': Vp})

        # Map to appropriate space
        factor = sym.exp(- self.β / sym.Rational(2) * (Vq + Vp))

        # Mapped operator
        backward = eq.map_operator(forward, f, factor)

        # Discretize factor
        factor = quad.discretize(factor)

        return quad, forward, backward, factor

    def solve(self, backward, quad, factor, degrees):

        # Discretization of the operator
        mat = quad.discretize_op(backward, self.degree).matrix

        solutions = []
        for d in degrees:
            sub_mat = (mat[0:d+1, 0:d+1])
            if type(sub_mat) is np.ndarray:
                sub_mat = sub_mat.copy(order='C')
            eig_vals, eig_vecs = las.eigs(sub_mat, k=1, which='LR')
            ground_state = np.real(eig_vecs.T[0])
            ground_state = ground_state * np.sign(ground_state[0])
            ground_state_eval = quad.eval(quad.series(ground_state))*factor
            norm = quad.norm(ground_state_eval, n=1, flat=True)
            ground_state_eval = ground_state_eval / norm
            solutions.append(ground_state_eval)
        return solutions

    def test_exact_sol_gaussian(self):

        Vp, m, s2 = self.x*self.x/4, sym.Rational(1, 10), sym.Rational(1, 5)
        quad, forward, backward, factor = self.sym_calc(Vp, m, s2)
        solution = eq.solve_gaussian(forward, self.f, [self.x])
        norm_sol = quad.integrate(solution, flat=True)
        self.assertTrue(abs(norm_sol - 1) < 1e-6)

    def test_gaussian(self):

        Vp, m, s2 = self.x*self.x/4, sym.Rational(1, 10), sym.Rational(1, 5)
        quad, forward, backward, factor = self.sym_calc(Vp, m, s2)

        # Exact Solution
        solution = eq.solve_gaussian(forward, self.f, [self.x])
        solution_eval = quad.discretize(solution)

        # Numerical solutions
        degrees = list(range(5, self.degree))
        solutions = self.solve(backward, quad, factor, degrees)

        # Associated errors
        errors = []
        for sol in solutions:
            error = quad.norm(sol - solution_eval, flat=True)
            errors.append(error)

        log_errors = np.log(errors)
        poly_approx = np.polyfit(degrees, log_errors, 1)
        errors_approx = np.exp(np.polyval(poly_approx, degrees))
        error = la.norm(log_errors - np.log(errors_approx), 2)

        plt.semilogy(degrees, errors, 'k.')
        plt.semilogy(degrees, errors_approx)
        # plt.show()

        self.assertTrue(errors[-1] < 1e-12)
        self.assertTrue(error < 5)

    def test_norm_solution_bistable(self):

        # Parameters of potential associated with Hermite polynomials
        Vp = self.x**4/4 - self.x**2/2
        m, s2 = sym.Rational(1, 10), sym.Rational(1, 10)
        quad, forward, backward, factor = self.sym_calc(Vp, m, s2)

        sol_an = sym.exp(-self.β * Vp)
        sol_num = quad.discretize(sol_an)
        norm_1_an = quad.norm(sol_an, n=1, flat=True)
        norm_1_num = quad.norm(sol_num, n=1, flat=True)
        norm_2_an = quad.norm(sol_an, n=2, flat=True)
        norm_2_num = quad.norm(sol_num, n=2, flat=True)
        self.assertAlmostEqual(norm_1_an, norm_1_num)
        self.assertAlmostEqual(norm_2_an, norm_2_num)

    def test_bistable(self):

        # Parameters of potential associated with Hermite polynomials
        Vp = self.x**4/4 - self.x**2/2
        m, s2 = sym.Rational(1, 10), sym.Rational(1, 10)
        quad, forward, backward, factor = self.sym_calc(Vp, m, s2)

        # Exact Solution
        exact_sol = sym.exp(-self.β * Vp)
        norm = quad.norm(exact_sol, n=1, flat=True)
        exact_sol = exact_sol / norm
        solution_eval = quad.discretize(exact_sol)

        degrees = list(range(10, self.degree))
        solutions = self.solve(backward, quad, factor, degrees)

        # Associated errors
        errors = []
        for sol in solutions:
            error = quad.norm(sol - solution_eval, flat=True)
            errors.append(error)

        log_errors = np.log(errors)
        poly_approx = np.polyfit(degrees, log_errors, 1)
        errors_approx = np.exp(np.polyval(poly_approx, degrees))
        error = la.norm(log_errors - np.log(errors_approx), 2)

        plt.semilogy(degrees, errors, 'k.')
        plt.semilogy(degrees, errors_approx)
        # plt.show()

        self.assertTrue(errors[-1] < 1e-10)
        self.assertTrue(error < 10)

    def test_gaussian_sparse(self):
        hm.settings['sparse'] = True
        self.test_gaussian()

    def test_bistable_sparse(self):
        hm.settings['sparse'] = True
        self.test_bistable()


class TestConvergenceFokkerPlanck2d(unittest.TestCase):

    settings = {
        'cache': False,
        'cachedir': '/tmp/test_hermite',
        'tensorize': False,
        'sparse': False,
        'trails': False,
        'debug': False,
        }

    def setUp(self):

        # Equation parameters
        self.eq = eq.McKean_Vlasov
        self.x, self.y, self.f = self.eq.x, self.eq.y, self.eq.f

        # Default degree
        self.degree = 50

        # Set default settings
        hm.settings.update(self.settings)

    def test_correctness_operator(self):
        params = self.eq.params()
        β, γ, ε, θ, m = (params[x] for x in ['β', 'γ', 'ε', 'θ', 'm'])
        d, x, y, f = sym.diff, self.x, self.y, self.f
        Vp, Vy = params['Vp'], params['Vy']
        explicit = d(d(Vp, x)*f + θ*(x-m)*f - (1-γ)*sym.sqrt(1/β)*y*f/ε, x) \
            + γ**2/β * d(d(f, x), x) \
            + (1/ε**2) * d(Vy.diff(y) * f, y) \
            + (1/ε**2) * d(d(f, y), y)
        self.assertTrue((explicit - self.eq.equation(params)).cancel() == 0)

    def sym_calc(self, Vp, parameters, m, s2x, s2y,
                 degree=None):

        # Defalut degree
        degree = degree or self.degree

        # Number of quadrature points
        n_points_num = 2*degree + 1

        # Equation parameters
        β = parameters['β']

        equation = eq.McKean_Vlasov
        x, y, f = equation.x, equation.y, equation.f

        # Calculation of the solution
        new_q = hm.Quad.gauss_hermite
        quad = new_q(n_points_num, dim=2,
                     mean=[m, 0], cov=[[s2x, 0], [0, s2y]])

        # Potential for approximation
        Vqx = sym.Rational(1/2)*(x-m)*(x-m)/(β*s2x)
        Vqy = sym.Rational(1/2)*y*y/s2y

        # Fokker Planck for McKean-Vlasov equation
        parameters.update({'Vp': Vp, 'Vy': y*y/2})
        forward = equation.equation(parameters)

        # Map to appropriate space
        factor_x = sym.exp(- β / 2 * (Vqx + Vp))
        factor_y = sym.exp(- 1/2 * (Vqy + y*y/2))
        factor = factor_x * factor_y

        # Mapped operator
        backward = eq.map_operator(forward, f, factor)

        # Discretize factor
        factor_x = quad.project(0).discretize(factor_x)
        factor_y = quad.project(1).discretize(factor_y)
        factor = quad.discretize(factor)

        return quad, forward, backward, factor, factor_x, factor_y

    def solve(self, backward, quad, factor, degrees):

        # Discretization of the operator
        mat = quad.discretize_op(backward, degrees[-1]).matrix

        solutions = []

        v0, eig_vec = None, None
        for d in degrees:
            # print(d)
            npolys = core.iterator_size(2, d)
            if d is not degrees[0]:
                v0 = np.zeros(npolys)
                for i in range(len(eig_vec)):
                    v0[i] = eig_vec[i]
            sub_mat = (mat[0:npolys, 0:npolys])
            if type(sub_mat) is np.ndarray:
                sub_mat = sub_mat.copy(order='C')
            # import ipdb; ipdb.set_trace()
            eig_vals, eig_vecs = las.eigs(sub_mat, k=1, v0=v0, which='LR')
            eig_vec = np.real(eig_vecs.T[0])
            ground_state = eig_vec * np.sign(eig_vec[0])
            ground_state_eval = quad.eval(quad.series(ground_state))*factor
            norm = quad.norm(ground_state_eval, n=1, flat=True)
            ground_state_eval = ground_state_eval / norm
            solutions.append(ground_state_eval)

            # fig, ax = plt.subplots(1, 1)
            # quad.plot(quad.series(ground_state), factor, ax=ax)
            # plt.show()

        return solutions, quad.series(ground_state)

    def test_exact_sol_gaussian(self):

        r = sym.Rational
        params = {'β': r(3), 'ε': r(1), 'γ': 0, 'θ': 0, 'm': 0}
        Vp = self.x*self.x/4
        quad, forward, _, factor, _, _ = self.sym_calc(Vp, params, 0, 1, 1)
        solution = eq.solve_gaussian(forward, self.f, [self.x, self.y])
        norm_sol = quad.norm(solution, n=1, flat=True)
        self.assertTrue(abs(norm_sol - 1) < 1e-6)

    def test_gaussian(self):

        r = sym.Rational
        Vp, m, s2x, s2y = self.x**2/4, r(1, 10), r(1, 2), 1/sym.sqrt(2)
        params = {'β': r(1), 'ε': r(1), 'γ': 0, 'θ': 0, 'm': 0}
        args = [Vp, params, m, s2x, s2y]
        quad, forward, backward, factor, _, _ = self.sym_calc(*args)

        # Exact Solution
        solution = eq.solve_gaussian(forward, self.f, [self.x, self.y])
        solution_eval = quad.discretize(solution)

        # Numerical solutions
        degrees = list(range(10, self.degree))
        solutions, finest = self.solve(backward, quad, factor, degrees)

        # Plot of the finest solution
        # fig, ax = plt.subplots(1, 1)
        # quad.plot(finest, factor, ax=ax)
        # plt.show()

        # Associated errors
        errors = []
        for sol in solutions:
            error = quad.norm(sol - solution_eval, flat=True)
            errors.append(error)
            # print(error)

        log_errors = np.log(errors)
        poly_approx = np.polyfit(degrees, log_errors, 1)
        errors_approx = np.exp(np.polyval(poly_approx, degrees))
        error = la.norm(log_errors - np.log(errors_approx), 2)

        plt.semilogy(degrees, errors, 'k.')
        plt.semilogy(degrees, errors_approx)
        # plt.show()

        self.assertTrue(errors[-1] < 1e-6)
        self.assertTrue(error < 1)

    def test_bistable(self):

        r = sym.Rational
        Vp, degree = self.x**4/4 - self.x**2/2, 40
        m, s2x, s2y = r(1, 10), r(1, 20), 1
        params = {'β': r(10), 'ε': r(.5), 'γ': 0, 'θ': 0, 'm': 0}
        args = [Vp, params, m, s2x, s2y, degree]
        quad, forward, backward, factor, _, _ = self.sym_calc(*args)

        # Numerical solutions
        degrees = list(range(5, degree))
        solutions, finest = self.solve(backward, quad, factor, degrees)
        finest_eval = solutions[-1]

        # Plot of the finest solution
        # fig, ax = plt.subplots(1, 1)
        # quad.plot(finest, factor, ax=ax)
        # plt.show()

        # Associated errors
        errors, degrees = [], degrees[0:-1]
        for sol in solutions[0:-1]:
            error = quad.norm(sol - finest_eval, flat=True)
            errors.append(error)
            # print(error)

        log_errors = np.log(errors)
        poly_approx = np.polyfit(degrees, log_errors, 1)
        errors_approx = np.exp(np.polyval(poly_approx, degrees))
        error = la.norm(log_errors - np.log(errors_approx), 2)

        # plt.semilogy(degrees, errors, 'k.')
        # plt.semilogy(degrees, errors_approx)
        # plt.show()

        self.assertTrue(errors[-1] < 1e-3)
        self.assertTrue(error < 3)

    def test_gaussian_sparse(self):
        hm.settings['sparse'] = True
        self.test_gaussian()

    def test_gaussian_tensorize(self):
        hm.settings['tensorize'] = True
        self.test_gaussian()

    def test_gaussian_sparse_tensorize(self):
        hm.settings['sparse'] = True
        hm.settings['tensorize'] = True
        self.test_gaussian()

    def test_bistable_sparse(self):
        hm.settings['sparse'] = True
        self.test_bistable()

    def test_bistable_tensorize(self):
        hm.settings['tensorize'] = True
        self.test_bistable()

    def test_bistable_sparse_tensorize(self):
        hm.settings['sparse'] = True
        hm.settings['tensorize'] = True
        self.test_bistable()


# class TestConvergenceFokkerPlanck3d(unittest.TestCase):

#     def setUp(self):

#         # Equation parameters
#         self.eq = eq.McKean_Vlasov_harmonic_noise
#         self.x, self.y, self.z = self.eq.x, self.eq.y, self.eq.z
#         self.f = self.eq.f
#         self.dim = 3

#     def test_correctness_operator(self):
#         params = self.eq.params()
#         d, x, y, z, f = sym.diff, self.x, self.y, self.z, self.f
#         β, γ, ε, θ, m = (params[x] for x in ['β', 'γ', 'ε', 'θ', 'm'])
#         λ = 1
#         Vp = params['Vp']
#         explicit = d(d(Vp, x)*f + θ*(x-m)*f - (1-γ)*sym.sqrt(1/β)*z*f/ε, x) \
#             + γ**2/β * d(d(f, x), x) \
#             + (1/ε**2) * (-y*d(f, z) + z*d(f, y)) \
#             + (1/ε**2) * λ*d(f*y + d(f, y), y)
#         self.assertTrue((explicit - self.eq.equation(params)).cancel() == 0)

#     def sym_calc(self, Vp, parameters, s2x, s2y, s2z, degree=10):

#         # Number of quadrature points
#         n_points_num = 2*degree + 1

#         # Equation parameters
#         β = parameters['β']

#         equation = eq.McKean_Vlasov_harmonic_noise
#         x, y, z, f = self.x, self.y, self.z, self.f

#         # Calculation of the solution
#         new_q = hm.Quad.gauss_hermite
#         cov = [[s2x, 0, 0], [0, s2y, 0], [0, 0, s2z]]
#         quad = new_q(n_points_num, dim=3, mean=[0]*3, cov=cov)

#         # Potential for approximation
#         Vqx = sym.Rational(1/2)*x*x/(β*s2x)
#         Vqy = sym.Rational(1/2)*y*y/s2y
#         Vqz = sym.Rational(1/2)*z*z/s2z

#         # Fokker Planck for McKean-Vlasov equation
#         parameters.update({'Vp': Vp})
#         forward = equation.equation(parameters)

#         # Map to appropriate space
#         factor_x = sym.exp(- β / 2 * (Vqx + Vp))
#         factor_y = sym.exp(- 1/2 * (y*y/2 + Vqy))
#         factor_z = sym.exp(- 1/2 * (z*z/2 + Vqz))
#         factor = factor_x * factor_y * factor_z
#         factors = [factor_x, factor_y, factor_z]

#         # Mapped operator
#         backward = eq.map_operator(forward, f, factor)

#         # Probability density fluxes
#         fluxes = eq.McKean_Vlasov_harmonic_noise.fluxes(parameters)
#         fluxes = [eq.map_operator(flux, f, factor) for flux in fluxes]

#         return quad, backward, fluxes, factor, factors

#     def check_consistency_varf(self, sparse=True):
#         r = sym.Rational
#         Vp, degree = self.x**4/4 - self.x**2/2, 10
#         s2x, s2y, s2z = r(1, 2), 1, 1
#         params = {'β': 5, 'ε': 0.5, 'γ': 0, 'θ': 0, 'm': 0}
#         args = [Vp, params, s2x, s2y, s2z, degree]
#         quad, backward, _, _, _ = self.sym_calc(*args)

#         varf_cross = quad.discretize_op(backward, degree,
#                                         sparse=sparse, index_set="cross")
#         varf_trian = quad.discretize_op(backward, degree + 2,
#                                         sparse=sparse, index_set="triangle")
#         self.assertTrue(varf_trian.to_cross(degree) == varf_cross)

#     def test_consistency_varf_dense(self):
#         hm.settings['tensorize'] = True
#         self.check_consistency_varf(sparse=False)

#     def test_consistency_varf_sparse(self):
#         hm.settings['tensorize'] = True
#         self.check_consistency_varf(sparse=True)

#     def solve(self, backward, quad, factors, degrees,
#               index_set="triangle", fluxes=None):

#         # Discretization of the operator
#         hm.settings['tensorize'] = True
#         hm.settings['cache'] = True
#         hm.settings['cachedir'] = '/tmp/test_fp_3d'
#         hm.settings['trails'] = True
#         hm.settings['debug'] = False

#         kwargs = {'sparse': True, 'index_set': index_set}
#         var = quad.discretize_op(backward, degrees[-1], **kwargs)
#         fd = [quad.discretize_op(flux, degrees[-1], **kwargs) for flux in fluxes]
#         mat = var.matrix

#         solutions = []
#         dirs_removed = [1]
#         factor_removed, weight_removed = 1, 1
#         for d in dirs_removed:
#             factor_removed *= factors[d]
#             weight_removed *= quad.position.weights()[d]
#         f_projection = factor_removed / weight_removed
#         f_projection = sym.Rational(1)
#         s_projection = quad.project(dirs_removed).\
#             transform(f_projection, degrees[-1], index_set=index_set)
#         dirs_visu = [i for i in range(3) if i not in dirs_removed]
#         factor_visu = 1
#         for d in dirs_visu:
#             factor_visu *= factors[d]

#         # Quadrature for vizualization
#         nv, bx, by, bz = 200, 3, 5, 5
#         quad_visu = hm.Quad.newton_cotes([nv, nv, nv], [bx, by, bz])
#         quad_visu = quad_visu.project(dirs_visu)

#         import ipdb
#         ipdb.set_trace()

#         v0, eig_vec = None, None
#         for d in degrees:
#             npolys = core.iterator_size(self.dim, d, index_set=index_set)
#             if d is not degrees[0]:
#                 v0 = np.zeros(npolys)
#                 for i in range(len(eig_vec)):
#                     v0[i] = eig_vec[i]
#             sub_mat = (mat[0:npolys, 0:npolys])

#             if type(sub_mat) is np.ndarray:
#                 sub_mat = sub_mat.copy(order='C')
#             eigs = stats.log_stats()(cache.cache(quiet=True)(las.eigs))
#             eig_vals, eig_vecs = eigs(sub_mat, k=1, v0=v0, which='LR')
#             eig_vec = np.real(eig_vecs.T[0])
#             ground_state = eig_vec * np.sign(eig_vec[0])
#             ground_state_series = quad.series(ground_state,
#                                               index_set=index_set)
#             # ground_state_eval = quad.eval(ground_state_series)*factor
#             # norm = quad.norm(ground_state_eval, n=1, flat=True)
#             # ground_state_eval = ground_state_eval / norm
#             # solutions.append(ground_state_eval)

#             fluxes_eq = [flux.subdegree(d)(ground_state_series) for flux in fd]
#             sub_projection = s_projection.subdegree(d)
#             inner_series = ground_state_series.inner(sub_projection)
#             flux_proj = [flux.inner(sub_projection) for flux in fluxes_eq]

#             fig, ax = plt.subplots(1, 2)
#             cont = quad_visu.plot(inner_series, factor_visu, ax=ax[0])
#             streams = quad_visu.\
#                 streamlines(flux_proj[0], flux_proj[1], factor_visu, ax=ax[0],
#                             colors='k', vmax=-.2, vmin=.2)
#             inner_series.plot(ax=ax[1])
#             plt.colorbar(cont, ax=ax[0])

#             fig, ax = plt.subplots(1, 2)
#             # factor_visu_d = quad_visu.discretize(factor_visu)
#             # fx, fy = quad_visu.eval(flux_proj[0]), quad_visu.eval(flux_proj[1])
#             cont1 = quad_visu.plot(flux_proj[0], factor_visu, ax=ax[0])
#             cont2 = quad_visu.plot(flux_proj[1], factor_visu, ax=ax[1])
#             plt.colorbar(cont1, ax=ax[0])
#             plt.colorbar(cont2, ax=ax[1])
#             plt.show()

#             # quad_visu.quiver(flux_proj[0], flux_proj[1],
#             #                  factor=factor_visu, ax=ax)
#             fig, ax = plt.subplots(1)
#             streams = quad_visu.streamlines(flux_proj[0], flux_proj[1],
#                                             factor_visu, ax)
#             plt.colorbar(streams, ax=ax)
#             plt.show()

#             import ipdb
#             ipdb.set_trace()

#         return solutions, quad.series(ground_state)

#     def test_bistable_cross(self):

#         r = sym.Rational
#         Vp, degree = self.x**4/4 - self.x**2/2, 80
#         s2x, s2y, s2z = r(1, 2), r(1, 1), r(1, 1)
#         params = {'β': r(1, 1), 'ε': r(1, 1), 'γ': 0, 'θ': r(2), 'm': 0}
#         args = [Vp, params, s2x, s2y, s2z, degree]
#         quad, forward, backward, factor, fx, fy, fz = self.sym_calc(*args)

#         # Numerical solutions
#         degrees = list(range(50, degree + 1, 5))
#         solutions, finest = self.solve(backward, quad, [fx, fy, fz], degrees,
#                                        index_set="cross")

#     def test_bistable(self):

#         r = sym.Rational
#         Vp, degree = self.x**4/4 - self.x**2/2, 40
#         s2x, s2y, s2z = r(1, 8), r(1, 1), r(1, 1)
#         params = {'β': r(1), 'ε': r(1, 2), 'γ': 0, 'θ': .5, 'm': .0}
#         args = [Vp, params, s2x, s2y, s2z, degree]
#         quad, backward, fluxes, factor, factors = self.sym_calc(*args)

#         # Numerical solutions
#         degrees = list(range(20, degree + 1, 5))
#         solutions, finest = self.solve(backward, quad, factors, degrees,
#                                        index_set="triangle", fluxes=fluxes)
#         finest_eval = solutions[-1]

#         # Plot of the finest solution
#         fig, ax = plt.subplots(1, 1)
#         quad.plot(finest, factor, ax=ax)
#         plt.show()

#         # Associated errors
#         errors, degrees = [], degrees[0:-1]
#         for sol in solutions[0:-1]:
#             error = quad.norm(sol - finest_eval, flat=True)
#             errors.append(error)
#             print(error)

#         log_errors = np.log(errors)
#         poly_approx = np.polyfit(degrees, log_errors, 1)
#         errors_approx = np.exp(np.polyval(poly_approx, degrees))
#         error = la.norm(log_errors - np.log(errors_approx), 2)

#         plt.semilogy(degrees, errors, 'k.')
#         plt.semilogy(degrees, errors_approx)
#         plt.show()

#         self.assertTrue(errors[-1] < 1e-3)
#         self.assertTrue(error < 1)


if __name__ == '__main__':
    unittest.main()
