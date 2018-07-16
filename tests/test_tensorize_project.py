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

import hermipy.core as core
import hermipy as hm

import unittest
import sympy as sym
import numpy.linalg as la
import scipy.sparse.linalg as las


class TestProject(unittest.TestCase):

    def setUp(self):
        hm.settings['tensorize'] = True
        hm.settings['cache'] = False

    def test_project_vector(self):
        n_points = 200
        degree = 10
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        coeffs_1d = quad_1d.transform('exp(x)', degree).coeffs
        coeffs_2d = quad_2d.transform('exp(x)', degree).coeffs
        projection = core.project(coeffs_2d, 2, 0)
        diff = la.norm(coeffs_1d - projection, 2)
        self.assertAlmostEqual(diff, 0)

    def test_project_matrix(self):
        n_points, degree = 200, 10
        function = 'exp(x)'
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        varf_1d = quad_1d.varf(function, degree).matrix
        varf_2d = quad_2d.varf(function, degree).matrix
        projection = core.project(varf_2d, 2, 0)
        diff = (la.norm(varf_1d - projection, 2))
        self.assertAlmostEqual(diff, 0)

    def test_project_matrix_2d_subspace(self):
        n_points, degree = 50, 10
        function = 'exp(x) * cos(y)'
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        quad_3d = hm.Quad.gauss_hermite(n_points, dim=3)
        varf_2d = quad_2d.varf(function, degree).matrix
        varf_3d = quad_3d.varf(function, degree).matrix
        projection = core.project(varf_3d, 3, [0, 1])
        diff = (la.norm(varf_2d - projection, 2))
        self.assertAlmostEqual(diff, 0)

    def test_project_matrix_2d_subspace_sparse(self):
        # hm.settings['tensorize'] = False
        n_points, degree = 50, 10
        function = 'exp(x) * cos(y)'
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        quad_3d = hm.Quad.gauss_hermite(n_points, dim=3)
        varf_2d = quad_2d.varf(function, degree, sparse=True).matrix
        varf_3d = quad_3d.varf(function, degree, sparse=True).matrix
        projection = core.project(varf_3d, 3, [0, 1])
        diff = (las.norm(varf_2d - projection))
        self.assertAlmostEqual(diff, 0)


class TestTensorize(unittest.TestCase):

    def setUp(self):
        hm.settings['tensorize'] = False
        hm.settings['cache'] = False

    def test_tensorize_vector(self):
        n_points = 200
        degree = 10
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        coeffs_1d = quad_1d.transform('exp(x)', degree).coeffs
        coeffs_2d = quad_2d.transform('exp(x)', degree).coeffs
        tensorized_coeffs_1d = core.tensorize(coeffs_1d, 2, 0)
        self.assertAlmostEqual(la.norm(coeffs_2d - tensorized_coeffs_1d, 2), 0)

    def test_tensorize_matrix(self):
        n_points = 200
        degree = 10
        function = 'exp(x)'
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        varf_1d = quad_1d.varf(function, degree).matrix
        varf_2d = quad_2d.varf(function, degree).matrix
        tensorized_varf_1d = core.tensorize(varf_1d, 2, 0)
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
        tensorized_coeffs = core.tensorize([coeffs_x, coeffs_y])
        diff = la.norm(coeffs_2d - tensorized_coeffs, 2)
        self.assertAlmostEqual(diff, 0)

    def test_tensorize_matrices(self):
        n_points = 200
        degree = 10
        f, fx, fy = 'exp(x) * (y*y*y)', 'exp(x)', 'x*x*x'
        quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        varf_x = quad_1d.varf(fx, degree).matrix
        varf_y = quad_1d.varf(fy, degree).matrix
        varf_2d = quad_2d.varf(f, degree).matrix
        tensorized_varf = core.tensorize([varf_x, varf_y])
        diff = la.norm(varf_2d - tensorized_varf, 2)
        self.assertAlmostEqual(diff, 0)

    def test_tensorize_quad(self):
        n_points = 200
        mean, cov = [.1, .2], [[1, 0], [0, 2]]
        quad = hm.Quad.gauss_hermite(n_points, dim=2, mean=mean, cov=cov)
        quad_x, quad_y = quad.project(0), quad.project(1)
        self.assertTrue(quad_x * quad_y == quad)

    def test_tensorize_series(self):
        n_points, degree = 200, 10
        mean, cov = [.1, .2], [[1, 0], [0, 2]]
        quad = hm.Quad.gauss_hermite(n_points, dim=2, mean=mean, cov=cov)
        series = quad.transform('exp(x) * (y*y*y)', degree)
        series_x, series_y = series.project(0), series.project(1)
        series_y = series_y * (1/series_y.coeffs[0])
        self.assertTrue(series_x * series_y == series)

    def test_tensorize_varf(self):
        n_points = 200
        degree = 10
        f, fx, fy = 'exp(x) * (y*y*y)', 'exp(x)', 'y*y*y'
        quad = hm.Quad.gauss_hermite(n_points, dim=2)
        varf_x = quad.project(0).varf(fx, degree)
        varf_y = quad.project(1).varf(fy, degree)
        varf_xy = quad.varf(f, degree)
        self.assertTrue(varf_x * varf_y == varf_xy)

    def test_tensorize_indices_yz(self):
        n_points, degree = 50, 10
        quad = hm.Quad.gauss_hermite(n_points, dim=3)
        quad_yz = quad.project([1, 2])
        y, z = hm.y, hm.z
        fyz = sym.cos(y)*z*z*z*z
        varf_yz = quad_yz.varf(fyz, degree)
        varf_yz_t = quad_yz.varf(fyz, degree, tensorize=True)
        self.assertTrue(varf_yz == varf_yz_t)

    def test_tensorize_matrix_2d_subspace(self):
        n_points, degree = 50, 10
        x, y, z = hm.x, hm.y, hm.z
        fx, fyz = sym.exp(x), sym.cos(y)*z*z*z*z
        f = fx * fyz
        quad = hm.Quad.gauss_hermite(n_points, dim=3)
        quad_x = quad.project([0])
        quad_yz = quad.project([1, 2])
        varf_x = quad_x.varf(fx, degree)
        varf_yz = quad_yz.varf(fyz, degree, tensorize=False)
        varf_xyz = quad.varf(f, degree, tensorize=True)
        tensorized = hm.Varf.tensorize([varf_x, varf_yz])
        self.assertTrue(varf_xyz == tensorized)
