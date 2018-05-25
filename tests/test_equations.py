import os
import pdb
import unittest
import sympy as sym
import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as las
import matplotlib.pyplot as plt

import hermite.quad as hm
import hermite.equations as eq
import hermite.settings as rc

from scipy.special import binom

settings = {
        'cache': True,
        'cachedir': '/tmp/test_hermite',
        'tensorize': False
        }
rc.settings.update(settings)
if not os.path.exists(settings['cachedir']):
    os.makedirs(settings['cachedir'])


class TestConvergenceFokkerPlanck1d(unittest.TestCase):

    def setUp(self):

        # Equation parameters
        equation = eq.Fokker_Planck_1d
        self.x, self.f = equation.x, equation.f

        # Inverse temperature
        self.β = sym.Rational(3)

        # Common parameters of numerical approximation
        self.degree = 100
        self.n_points_num = 2*self.degree + 1

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
        mat = quad.discretize_op(backward, self.f, self.degree, 2)

        solutions = []
        for d in degrees:
            sub_mat = (mat[0:d+1, 0:d+1]).copy(order='C')
            eig_vals, eig_vecs = las.eigs(sub_mat, k=1, which='LR')
            ground_state = np.real(eig_vecs.T[0])
            ground_state = ground_state * np.sign(ground_state[0])
            ground_state_eval = quad.eval(quad.series(ground_state))*factor
            norm = quad.norm(ground_state_eval, n=1, l2=True)
            ground_state_eval = ground_state_eval / norm
            solutions.append(ground_state_eval)
        return solutions

    def test_exact_sol_gaussian(self):

        Vp, m, s2 = self.x*self.x/4, sym.Rational(1, 10), sym.Rational(1, 5)
        quad, forward, backward, factor = self.sym_calc(Vp, m, s2)
        solution = eq.solve_gaussian(forward, self.f, [self.x])
        norm_sol = quad.integrate(solution, l2=True)
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
            error = quad.norm(sol - solution_eval, l2=True)
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
        norm_1_an = quad.norm(sol_an, n=1, l2=True)
        norm_1_num = quad.norm(sol_num, n=1, l2=True)
        norm_2_an = quad.norm(sol_an, n=2, l2=True)
        norm_2_num = quad.norm(sol_num, n=2, l2=True)
        self.assertAlmostEqual(norm_1_an, norm_1_num)
        self.assertAlmostEqual(norm_2_an, norm_2_num)

    def test_bistable(self):

        # Parameters of potential associated with Hermite polynomials
        Vp = self.x**4/4 - self.x**2/2
        m, s2 = sym.Rational(1, 10), sym.Rational(1, 10)
        quad, forward, backward, factor = self.sym_calc(Vp, m, s2)

        # Exact Solution
        exact_sol = sym.exp(-self.β * Vp)
        norm = quad.norm(exact_sol, n=1, l2=True)
        exact_sol = exact_sol / norm
        solution_eval = quad.discretize(exact_sol)

        degrees = list(range(10, self.degree))
        solutions = self.solve(backward, quad, factor, degrees)

        # Associated errors
        errors = []
        for sol in solutions:
            error = quad.norm(sol - solution_eval, l2=True)
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


class TestConvergenceFokkerPlanck2d(unittest.TestCase):

    def setUp(self):

        # Equation parameters
        equation = eq.McKean_Vlasov
        self.x, self.y, self.f = equation.x, equation.y, equation.f

    def sym_calc(self, Vp, parameters, m, s2x, s2y,
                 degree=50):

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
        parameters.update({'Vp': Vp})
        forward = equation.equation(parameters)

        # Map to appropriate space
        factor_x = sym.exp(- β / 2 * (Vqx + Vp))
        factor_y = sym.exp(- 1/2 * (Vqy + sym.sqrt(2)*y*y/2))
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
        mat = quad.discretize_op(backward, self.f, degrees[-1], 2)

        solutions = []

        v0, eig_vec = None, None
        for d in degrees:
            print(d)
            npolys = int(binom(d + 2, d))
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
            norm = quad.norm(ground_state_eval, n=1, l2=True)
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
        norm_sol = quad.norm(solution, n=1, l2=True)
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
        fig, ax = plt.subplots(1, 1)
        quad.plot(finest, factor, ax=ax)
        # plt.show()

        # Associated errors
        errors = []
        for sol in solutions:
            error = quad.norm(sol - solution_eval, l2=True)
            errors.append(error)
            # print(error)

        log_errors = np.log(errors)
        poly_approx = np.polyfit(degrees, log_errors, 1)
        errors_approx = np.exp(np.polyval(poly_approx, degrees))
        error = la.norm(log_errors - np.log(errors_approx), 2)

        plt.semilogy(degrees[0:-1], errors, 'k.')
        plt.semilogy(degrees[0:-1], errors_approx)
        # plt.show()

        self.assertTrue(errors[-1] < 1e-6)
        self.assertTrue(error < 1)

    def test_bistable(self):

        r = sym.Rational
        Vp, degree = self.x**4/4 - self.x**2/2, 40
        m, s2x, s2y = r(1, 10), r(1, 20), 1/2
        params = {'β': r(10), 'ε': r(.5), 'γ': 0, 'θ': 0, 'm': 0}
        args = [Vp, params, m, s2x, s2y, degree]
        quad, forward, backward, factor, _, _ = self.sym_calc(*args)

        # Numerical solutions
        degrees = list(range(5, degree))
        solutions, finest = self.solve(backward, quad, factor, degrees)
        finest_eval = solutions[-1]

        # Plot of the finest solution
        fig, ax = plt.subplots(1, 1)
        quad.plot(finest, factor, ax=ax)
        plt.show()

        # Associated errors
        errors, degrees = [], degrees[0:-1]
        for sol in solutions[0:-1]:
            error = quad.norm(sol - finest_eval, l2=True)
            errors.append(error)
            print(error)

        log_errors = np.log(errors)
        poly_approx = np.polyfit(degrees, log_errors, 1)
        errors_approx = np.exp(np.polyval(poly_approx, degrees))
        error = la.norm(log_errors - np.log(errors_approx), 2)

        plt.semilogy(degrees, errors, 'k.')
        plt.semilogy(degrees, errors_approx)
        plt.show()

        self.assertTrue(errors[-1] < 1e-3)
        self.assertTrue(error < 1)
