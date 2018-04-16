from . import hermite as hm
import unittest
import numpy as np
import numpy.polynomial.hermite_e as herm
import numpy.linalg as la
import math


class TestIntegrate(unittest.TestCase):

    def test_normalization_nodes(self):
        deg = [2**i for i in range(8)]
        for i in deg:
            quad = hm.Quad.gauss_hermite(n_points=i, dim=1)
            self.assertAlmostEqual(quad.integrate('1'), 1)
            self.assertAlmostEqual(quad.integrate('v[0]'), 0)

    def test_normalization_dim(self):
        for i in range(1, 4):
            quad = hm.Quad.gauss_hermite(n_points=100, dim=i)
            self.assertAlmostEqual(quad.integrate('1'), 1)
            self.assertAlmostEqual(quad.integrate('v[0]'), 0)
            self.assertAlmostEqual(quad.integrate('v[0]*v[0]'), 1)

    def test_mean(self):
        dim = 3
        mean = np.random.random(dim)
        quad = hm.Quad.gauss_hermite(n_points=8, dim=dim, mean=mean)
        for i in range(dim):
            fun = 'v[{}]'.format(i)
            coord = quad.integrate(fun)
            self.assertAlmostEqual(coord, mean[i])

    def test_covariance(self):
        dim = 3
        rand_mat = np.random.random((dim, dim))
        cov = np.matmul(rand_mat.T, rand_mat)
        quad = hm.Quad.gauss_hermite(n_points=8, dim=dim, cov=cov)
        for i in range(len(cov)):
            for j in range(len(cov)):
                fun = 'v[{}]*v[{}]'.format(i, j)
                cov_ij = quad.integrate(fun)
                self.assertAlmostEqual(cov_ij, cov[i][j])

    def test_all(self):
        dim = 3
        rand_mat = np.random.random((dim, dim))
        mean = np.random.random(dim)
        cov = np.matmul(rand_mat.T, rand_mat)
        quad = hm.Quad.gauss_hermite(n_points=8, dim=dim, mean=mean, cov=cov)
        for i in range(len(cov)):
            mean_i = quad.integrate('v[{}]'.format(i))
            self.assertAlmostEqual(mean_i, mean[i])
            for j in range(len(cov)):
                fun = '(v[{}]-{})*(v[{}]-{})'.format(i, mean[i], j, mean[j])
                cov_ij = quad.integrate(fun)
                self.assertAlmostEqual(cov_ij, cov[i][j])


class TestHermiteTransform(unittest.TestCase):

    def test_constant(self):
        degree = 30
        quad = hm.Quad.gauss_hermite(n_points=[degree, degree, degree])
        coeffs = quad.transform('1', degree).coeffs
        for i in range(len(coeffs)):
            target_value = 1. if i == 0 else 0.
            self.assertAlmostEqual(coeffs[i], target_value)

    def test_consistent_eval(self):
        degree = 10
        n_points = 10
        quad = hm.Quad.gauss_hermite(n_points)
        nodes_scipy, weights_scipy = herm.hermegauss(n_points)
        for i in range(degree):
            coeffs = np.zeros(degree + 1)
            coeffs[i] = 1
            factor = math.sqrt(math.factorial(i))
            hi_nodes_scipy = herm.hermeval(nodes_scipy, coeffs)
            hi_nodes = quad.eval(coeffs, degree)
            diff = sum(abs(hi_nodes_scipy/factor - hi_nodes))
            self.assertAlmostEqual(diff, 0)

    def test_forward_backward(self):
        degree = 10
        quad = hm.Quad.gauss_hermite(degree + 1)
        f_hermite = np.random.random(degree + 1)
        f_grid = quad.eval(f_hermite, degree)
        f_hermite_new = quad.transform(f_grid, degree).coeffs
        diff = sum(abs(f_hermite - f_hermite_new))
        self.assertAlmostEqual(diff, 0)

    def test_eval_different_grids(self):
        n_points = 100
        degree = 50
        quad_1 = hm.Quad.gauss_hermite(n_points, mean=[2.], cov=[[2.]])
        quad_2 = hm.Quad.gauss_hermite(n_points, mean=[-1.], cov=[[.1]])
        series = quad_1.transform('x', degree)
        evaluation = quad_2.eval(series, degree)
        discretization = quad_2.discretize('x')
        self.assertAlmostEqual(la.norm(evaluation - discretization, 2), 0)

    def test_newton_cotes(self):
        mean = [2.]
        cov = [[.5]]
        degree = 10
        quad_1 = hm.Quad.gauss_hermite([200], mean=mean, cov=cov)
        quad_2 = hm.Quad.newton_cotes([10000], [6.], mean=mean, cov=cov)
        series1 = quad_1.transform('x', degree)
        series2 = quad_2.transform('x', degree)
        assert(la.norm(series1.coeffs - series2.coeffs, 2) < 1e-3)


class TestHermiteVarf(unittest.TestCase):

    def test_simple_varf(self):
        n_points = 100
        degree = 30
        quad = hm.Quad.gauss_hermite(n_points, dim=2)
        var = quad.varf('1', degree)
        self.assertAlmostEqual(la.norm(var - np.eye(len(var)), 2), 0)

    def test_simple_dvarf(self):
        n_points = 100
        degree = 10
        quad = hm.Quad.gauss_hermite(n_points)
        bk_ou = quad.varfd('1', degree, [0, 0]) - quad.varfd('x', degree, [0])
        off_diag = bk_ou - np.diag(np.diag(bk_ou))
        self.assertAlmostEqual(la.norm(off_diag, 2), 0)


class TestTensorize(unittest.TestCase):

    def test_vector_tensorize(self):
        n_points = 200
        degree = 10
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        coeffs_1d = quad_1d.transform('exp(x)', degree).coeffs
        coeffs_2d = quad_2d.transform('exp(x)', degree).coeffs
        tensorized_coeffs_1d = hm.tensorize(coeffs_1d, 2, 0)
        self.assertAlmostEqual(la.norm(coeffs_2d - tensorized_coeffs_1d, 2), 0)

    def test_matrix_tensorize(self):
        n_points = 200
        degree = 10
        function = 'exp(x)'
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        varf_1d = quad_1d.varf(function, degree)
        varf_2d = quad_2d.varf(function, degree)
        tensorized_varf_1d = hm.tensorize(varf_1d, 2, 0)
        diff = (la.norm(varf_2d - tensorized_varf_1d, 2))
        self.assertAlmostEqual(diff, 0)
