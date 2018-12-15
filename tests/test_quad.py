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

import hermipy as hm

import unittest
import numpy as np
import numpy.polynomial.hermite_e as herm
import numpy.linalg as la
import sympy as sym
import math
import os
import tempfile

import scipy.sparse as sp
import scipy.sparse.linalg as las

settings = {'cache': False, 'cachedir': '/tmp/test_hermite'}
hm.settings.update(settings)


class TestIntegrate(unittest.TestCase):

    def test_normalization_nodes(self):
        deg = [2**i for i in range(8)]
        for i in deg:
            quad = hm.Quad.gauss_hermite(n_points=i, dim=1)
            self.assertAlmostEqual(quad.integrate('1'), 1)
            self.assertAlmostEqual(quad.integrate('x0'), 0)

    def test_normalization_dim(self):
        for i in range(1, 4):
            quad = hm.Quad.gauss_hermite(n_points=100, dim=i)
            self.assertAlmostEqual(quad.integrate('1'), 1)
            self.assertAlmostEqual(quad.integrate('x0'), 0)
            self.assertAlmostEqual(quad.integrate('x0*x0'), 1)

    def test_mean(self):
        dim = 3
        mean = np.random.random(dim)
        quad = hm.Quad.gauss_hermite(n_points=8, dim=dim, mean=mean)
        for i in range(dim):
            fun = 'x{}'.format(i)
            coord = quad.integrate(fun)
            self.assertAlmostEqual(coord, mean[i])

    def test_covariance(self):
        dim = 3
        rand_mat = np.random.random((dim, dim))
        cov = np.matmul(rand_mat.T, rand_mat)
        quad = hm.Quad.gauss_hermite(n_points=8, dim=dim, cov=cov)
        for i in range(dim):
            for j in range(dim):
                fun = 'x{}*x{}'.format(i, j)
                cov_ij = quad.integrate(fun)
                self.assertAlmostEqual(cov_ij, cov[i][j])

    def test_all(self):
        dim = 3
        rand_mat = np.random.random((dim, dim))
        mean = np.random.random(dim)
        cov = np.matmul(rand_mat.T, rand_mat)
        quad = hm.Quad.gauss_hermite(n_points=8, dim=dim, mean=mean, cov=cov)
        for i in range(dim):
            mean_i = quad.integrate('x{}'.format(i))
            self.assertAlmostEqual(mean_i, mean[i])
            for j in range(dim):
                fun = '(x{}-{})*(x{}-{})'.format(i, mean[i], j, mean[j])
                cov_ij = quad.integrate(fun)
                self.assertAlmostEqual(cov_ij, cov[i][j])

    def test_l2_aligned(self):
        dim = 1
        rand_mat = np.random.random((dim, dim))
        mean = np.random.random(dim)
        cov = np.matmul(rand_mat.T, rand_mat)
        quad_aligned = hm.Quad.gauss_hermite(8, dim=dim, mean=mean, cov=cov)
        gaussian = quad_aligned.position.weight()
        integral = quad_aligned.integrate(gaussian, flat=True)
        self.assertAlmostEqual(integral, 1.)

    def test_l2_non_aligned(self):
        dim = 2
        cov = [[3, 1], [1, 6]]
        mean = [0.5, 1]
        quad_ali = hm.Quad.gauss_hermite(50, dim=dim, mean=mean, cov=cov)
        quad_std = hm.Quad.gauss_hermite(50, dim=dim)
        gaussian_ali = quad_ali.position.weight()
        gaussian_std = quad_std.position.weight()
        integral_1 = quad_ali.integrate(gaussian_std, flat=True)
        integral_2 = quad_std.integrate(gaussian_ali, flat=True)
        self.assertTrue(abs(integral_1 - 1.) < .01)
        self.assertTrue(abs(integral_2 - 1.) < .01)


class TestHermiteTransform(unittest.TestCase):

    def test_constant(self):
        degree = 30
        quad = hm.Quad.gauss_hermite(n_points=[degree, degree, degree])
        coeffs = quad.transform('1', degree).coeffs
        for i, c in enumerate(coeffs):
            target_value = 1. if i == 0 else 0.
            self.assertAlmostEqual(c, target_value)

    def test_consistent_eval(self):
        degree = 10
        n_points = 10
        quad = hm.Quad.gauss_hermite(n_points)
        nodes_scipy, _ = herm.hermegauss(n_points)
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
        self.assertTrue(la.norm(series1.coeffs - series2.coeffs, 2) < 1e-3)


class TestHermiteVarf(unittest.TestCase):

    def test_simple_varf(self):
        n_points = 100
        degree = 30
        quad = hm.Quad.gauss_hermite(n_points, dim=2)
        var = quad.varf('1', degree).matrix
        self.assertAlmostEqual(la.norm(var - np.eye(len(var)), 2), 0)

    def test_simple_dvarf(self):
        n_points = 100
        degree = 10
        quad = hm.Quad.gauss_hermite(n_points)
        bk_ou = quad.varfd('1', degree, [0, 0]).matrix \
            + quad.varfd('-x', degree, [0]).matrix
        off_diag = bk_ou - np.diag(np.diag(bk_ou))
        self.assertAlmostEqual(la.norm(off_diag, 2), 0)

    def test_varf_split(self):
        n_points = 100
        degree = 10
        quad = hm.Quad.gauss_hermite(n_points, dim=2)
        x, y = sym.symbols('x y')
        function = x*x*sym.cos(x) + sym.exp(y)*x + sym.sqrt(2) + 2
        v1 = quad.varf(function, degree, tensorize=False).matrix
        v2 = quad.varf(function, degree, tensorize=True).matrix
        diff = la.norm(v1 - v2, 2)
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
        dim = 5
        diag = 1 + np.abs(np.random.random(dim))
        self.cov = np.diag(diag)
        self.quad = hm.Quad.gauss_hermite(100, dim=dim, cov=self.cov)
        self.quad_low = hm.Quad.gauss_hermite(6, dim=dim, cov=self.cov)
        self.x = [sym.symbols('x' + str(i)) for i in range(dim)]
        hm.settings['tensorize'] = True
        hm.settings['trails'] = False

    def tearDown(self):
        hm.settings.update(settings)

    def test_integrate(self):
        dim = len(self.cov)
        for i in range(dim):
            for j in range(dim):
                fun = self.x[i]*self.x[j]
                cov_ij = self.quad.integrate(fun)
                self.assertAlmostEqual(cov_ij, self.cov[i][j])

    def test_sparse_varf_5d_simple(self):
        degree = 12
        function = 1
        sp_var = self.quad.varf(function, degree, sparse=True).matrix
        self.assertTrue(sp_var.nnz == sp_var.shape[0])

    def test_sparse_varf_5d(self):
        degree = 10
        function = 1. + self.x[0] + self.x[1] + self.x[1]*self.x[0] \
            + self.x[2]**2 + self.x[3]**3 + self.x[4]
        sp_var = self.quad.varf(function, degree, sparse=True).matrix
        self.assertTrue(sp_var.nnz < sp_var.shape[0] * sp_var.shape[1])

    def test_consistency(self):
        degree = 3
        function = 1. + self.x[0] + self.x[1] + self.x[1]*self.x[0] \
            + self.x[2]**2 + self.x[3]**3 + self.x[4]
        var_1 = self.quad_low.varf(function, degree, tensorize=True).matrix
        var_2 = self.quad_low.varf(function, degree, tensorize=False).matrix
        self.assertAlmostEqual(la.norm(var_1 - var_2), 0.)


class TestSparseFunctions(unittest.TestCase):

    def setUp(self):
        n_points = 100
        self.quad1 = hm.Quad.gauss_hermite(n_points)
        self.quad2 = hm.Quad.gauss_hermite(n_points, dim=2)
        self.x = [sym.symbols('x' + str(i)) for i in range(2)]
        hm.settings['tensorize'] = False
        hm.settings['cache'] = False

    def tearDown(self):
        hm.settings.update(settings)

    def test_sparse_varf_simple(self):
        degree = 30
        function = '1'
        sp_var = self.quad2.varf(function, degree, sparse=True).matrix
        self.assertEqual(sp_var.nnz, sp_var.shape[0])

    def test_sparse_varf1d(self):
        degree, degreef_max = 50, 7
        x = sym.symbols('x')
        function = 0
        for deg in range(degreef_max):
            function += x**deg
            sp_var = self.quad1.varf(function, degree, sparse=True).matrix
            coords = sp_var.tocoo()
            bw = 0
            for i, j, v in zip(coords.row, coords.col, coords.data):
                if abs(i - j) > bw:
                    bw = abs(i - j)
                    if bw > deg:
                        print(i, j, v)
                    self.assertTrue(bw <= deg)
            self.assertEqual(bw, deg)

    def test_sparse_varf_2d(self):
        degree = 50
        function = 1. + self.x[0] + self.x[1] + self.x[1]*self.x[0]
        sp_var = self.quad2.varf(function, degree, sparse=True)
        var = self.quad2.varf(function, degree, sparse=False)
        self.assertAlmostEqual(la.norm(var.matrix - sp_var.matrix, 2), 0)

    def test_sparse_varfd_1d(self):
        n_points = 100
        degree = 10
        quad = hm.Quad.gauss_hermite(n_points)
        mat1 = quad.varfd('1', degree, [0, 0], sparse=True).matrix
        mat2 = quad.varfd('x', degree, [0], sparse=True).matrix
        bk_ou = mat1 - mat2
        off_diag = bk_ou - sp.diags(bk_ou.diagonal())
        self.assertAlmostEqual(las.norm(off_diag), 0)


if __name__ == '__main__':
    unittest.main()
