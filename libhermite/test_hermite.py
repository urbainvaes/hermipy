from . import hermite as hm
import unittest
import numpy as np
import numpy.polynomial.hermite_e as herm
import numpy.linalg as la
import sympy as sym
import math
import os
import tempfile

settings = {'cache': False, 'cachedir': '/tmp/test_hermite'}
hm.settings.update(settings)
if not os.path.exists(settings['cachedir']):
    os.makedirs(settings['cachedir'])


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

    def test_l2_aligned(self):
        dim = 1
        rand_mat = np.random.random((dim, dim))
        mean = np.random.random(dim)
        cov = np.matmul(rand_mat.T, rand_mat)
        quad_aligned = hm.Quad.gauss_hermite(8, dim=dim, mean=mean, cov=cov)
        gaussian = quad_aligned.weight()
        integral = quad_aligned.integrate(gaussian, l2=True)
        self.assertAlmostEqual(integral, 1.)

    def test_l2_non_aligned(self):
        dim = 2
        cov = [[3, 1], [1, 6]]
        mean = [0.5, 1]
        quad_ali = hm.Quad.gauss_hermite(50, dim=dim, mean=mean, cov=cov)
        quad_std = hm.Quad.gauss_hermite(50, dim=dim)
        gaussian_ali = quad_ali.weight()
        gaussian_std = quad_std.weight()
        integral_1 = quad_ali.integrate(gaussian_std, l2=True)
        integral_2 = quad_std.integrate(gaussian_ali, l2=True)
        self.assertTrue(abs(integral_1 - 1.) < .01)
        self.assertTrue(abs(integral_2 - 1.) < .01)


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
            hi_nodes = quad.eval(coeffs)
            diff = sum(abs(hi_nodes_scipy/factor - hi_nodes))
            self.assertAlmostEqual(diff, 0)

    def test_forward_backward(self):
        degree = 10
        quad = hm.Quad.gauss_hermite(degree + 1)
        f_hermite = np.random.random(degree + 1)
        f_grid = quad.eval(f_hermite)
        f_hermite_new = quad.transform(f_grid, degree).coeffs
        diff = sum(abs(f_hermite - f_hermite_new))
        self.assertAlmostEqual(diff, 0)

    def test_eval_different_grids(self):
        n_points = 100
        degree = 50
        quad_1 = hm.Quad.gauss_hermite(n_points, mean=[2.], cov=[[2.]])
        quad_2 = hm.Quad.gauss_hermite(n_points, mean=[-1.], cov=[[.1]])
        series = quad_1.transform('x', degree)
        evaluation = quad_2.eval(series)
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

    def test_varf_split(self):
        n_points = 100
        degree = 10
        quad = hm.Quad.gauss_hermite(n_points, dim=2)
        x, y = sym.symbols('x y')
        function = x*x*sym.cos(x) + sym.exp(y)*x + sym.sqrt(2) + 2
        v1 = quad.varf(function, degree, tensorize=False)
        v2 = quad.varf(function, degree, tensorize=True)
        diff = la.norm(v1 - v2, 2)
        self.assertAlmostEqual(diff, 0)


class TestTensorize(unittest.TestCase):

    def test_tensorize_vector(self):
        n_points = 200
        degree = 10
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        coeffs_1d = quad_1d.transform('exp(x)', degree).coeffs
        coeffs_2d = quad_2d.transform('exp(x)', degree).coeffs
        tensorized_coeffs_1d = hm.tensorize(coeffs_1d, 2, 0)
        self.assertAlmostEqual(la.norm(coeffs_2d - tensorized_coeffs_1d, 2), 0)

    def test_tensorize_matrix(self):
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

    def test_tensorize_vectors(self):
        n_points = 200
        degree = 10
        f, fx, fy = 'exp(x) * (y*y*y)', 'exp(x)', 'x*x*x'
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        coeffs_x = quad_1d.transform(fx, degree).coeffs
        coeffs_y = quad_1d.transform(fy, degree).coeffs
        coeffs_2d = quad_2d.transform(f, degree).coeffs
        tensorized_coeffs = hm.tensorize([coeffs_x, coeffs_y])
        diff = la.norm(coeffs_2d - tensorized_coeffs, 2)
        self.assertAlmostEqual(diff, 0)

    def test_tensorize_matrices(self):
        n_points = 200
        degree = 10
        f, fx, fy = 'exp(x) * (y*y*y)', 'exp(x)', 'x*x*x'
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        varf_x = quad_1d.varf(fx, degree)
        varf_y = quad_1d.varf(fy, degree)
        varf_2d = quad_2d.varf(f, degree)
        tensorized_varf = hm.tensorize([varf_x, varf_y])
        diff = la.norm(varf_2d - tensorized_varf, 2)
        self.assertAlmostEqual(diff, 0)

    def test_project_vector(self):
        n_points = 200
        degree = 10
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        coeffs_1d = quad_1d.transform('exp(x)', degree).coeffs
        coeffs_2d = quad_2d.transform('exp(x)', degree).coeffs
        projection = hm.project(coeffs_2d, 2, 0)
        diff = la.norm(coeffs_1d - projection, 2)
        self.assertAlmostEqual(diff, 0)

    def test_project_matrix(self):
        n_points = 200
        degree = 10
        function = 'exp(x)'
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        varf_1d = quad_1d.varf(function, degree)
        varf_2d = quad_2d.varf(function, degree)
        projection = hm.project(varf_2d, 2, 0)
        diff = (la.norm(varf_1d - projection, 2))
        self.assertAlmostEqual(diff, 0)


class TestCache(unittest.TestCase):

    def setUp(self):
        x, y = sym.symbols('x y')
        n_points, dim = 100, 2
        self.quad = hm.Quad.gauss_hermite(n_points, dim=dim)
        self.function = x*x*sym.cos(x) + sym.exp(y)*x + sym.sqrt(2) + 2
        self.cachedir = tempfile.TemporaryDirectory()
        hm.settings['cachedir'] = self.cachedir.name
        hm.settings['cache'] = True

    def tearDown(self):
        self.cachedir.cleanup()
        hm.settings.update(settings)

    def test_varf(self):
        degree = 30
        self.quad.varf(self.function, degree)
        n_files_1 = len(os.listdir(self.cachedir.name))
        self.quad.varf(self.function, degree)
        n_files_2 = len(os.listdir(self.cachedir.name))
        self.assertEqual(n_files_1, n_files_2)

    def test_transform(self):
        degree = 30
        self.quad.transform(self.function, degree)
        n_files_1 = len(os.listdir(self.cachedir.name))
        self.quad.transform(self.function, degree)
        n_files_2 = len(os.listdir(self.cachedir.name))
        self.assertEqual(n_files_1, n_files_2)

    def test_integrate(self):
        self.quad.integrate(self.function)
        n_files_1 = len(os.listdir(self.cachedir.name))
        self.quad.integrate(self.function)
        n_files_2 = len(os.listdir(self.cachedir.name))
        self.assertEqual(n_files_1, n_files_2)


class TestTensorizeDecorator(unittest.TestCase):

    def setUp(self):
        n_points, dim = 100, 5
        diag = 1 + np.abs(np.random.random(dim))
        self.cov = np.diag(diag)
        self.quad = hm.Quad.gauss_hermite(n_points, dim=dim, cov=self.cov)
        self.v = [sym.symbols('v' + str(i)) for i in range(5)]
        hm.settings['tensorize'] = True

    def tearDown(self):
        hm.settings.update(settings)

    def test_integrate(self):
        for i in range(len(self.cov)):
            for j in range(len(self.cov)):
                fun = self.v[i]*self.v[j]
                cov_ij = self.quad.integrate(fun)
                self.assertAlmostEqual(cov_ij, self.cov[i][j])
