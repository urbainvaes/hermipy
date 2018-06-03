import hermite.core as core
import hermite.quad as hm
import unittest
import numpy.linalg as la
import scipy.sparse.linalg as las

import ipdb

class TestTensorize(unittest.TestCase):

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
        varf_1d = quad_1d.varf(function, degree)
        varf_2d = quad_2d.varf(function, degree)
        projection = core.project(varf_2d, 2, 0)
        diff = (la.norm(varf_1d - projection, 2))
        self.assertAlmostEqual(diff, 0)

    def test_project_matrix_2d_subspace(self):
        n_points, degree = 50, 10
        function = 'exp(x) * cos(y)'
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        quad_3d = hm.Quad.gauss_hermite(n_points, dim=3)
        varf_2d = quad_2d.varf(function, degree)
        varf_3d = quad_3d.varf(function, degree)
        projection = core.project(varf_3d, 3, ['x', 'y'])
        diff = (la.norm(varf_2d - projection, 2))
        self.assertAlmostEqual(diff, 0)

    def test_project_matrix_2d_subspace_sparse(self):
        n_points, degree = 50, 10
        function = 'exp(x) * cos(y)'
        quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
        quad_3d = hm.Quad.gauss_hermite(n_points, dim=3)
        varf_2d = quad_2d.varf(function, degree, sparse=True)
        varf_3d = quad_3d.varf(function, degree, sparse=True)
        projection = core.project(varf_3d, 3, ['x', 'y'])
        diff = (las.norm(varf_2d - projection))
        self.assertAlmostEqual(diff, 0)

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
        varf_1d = quad_1d.varf(function, degree)
        varf_2d = quad_2d.varf(function, degree)
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
        varf_x = quad_1d.varf(fx, degree)
        varf_y = quad_1d.varf(fy, degree)
        varf_2d = quad_2d.varf(f, degree)
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

#     def test_tensorize_varf(self):
#         n_points = 200
#         degree = 10
#         f, fx, fy = 'exp(x) * (y*y*y)', 'exp(x)', 'x*x*x'
#         quad = hm.Quad.gauss_hermite(n_points, dim=2)
#         varf_x = quad.project(0).varf(fx, degree)
#         varf_y = quad.project(1).varf(fy, degree)
#         varf_xy = quad.varf(f, degree)
#         self.assertTrue(varf_x * varf_y == varf_xy)

#         n_points, degree = 200, 10
#         mean, cov = [.1, .2], [[1, 0], [0, 2]]
#         quad = hm.Quad.gauss_hermite(n_points, dim=2, mean=mean, cov=cov)
#         varf = quad.transform('exp(x) * (y*y*y)', degree)
#         series_x, series_y = series.project(0), series.project(1)
#         series_y = series_y * (1/series_y.coeffs[0])
#         self.assertTrue(series_x * series_y == series)


#     # def test_tenorize_matrix_2d_subspace(self):
#     #     n_points, degree = 50, 10
#     #     x, y, z = sym.symbols('x y z', real=True)
#     #     fx, fyz = sym.exp(x), sym.cos(x)*y*y*y*y
#     #     f = fx * fyz.subs(y, z).subs(x, y)
#     #     quad = hm.Quad.gauss_hermite(n_points, dim=3)
#     #     quad_x = quad.project([0])
#     #     quad_xy = quad.project([1, 2])
#     #     varf_x = quad_x.varf(fx, degree)
#     #     varf_yz = quad_xy.varf(fyz, degree)
#     #     varf_xyz = quad.varf(f, degree)
#     #     # ipdb.set_trace()
#     #     tensorized = core.tensorize([varf_x, varf_yz], [[0], [1, 2]])
#     #     # diff = (la.norm(varf_2d - projection, 2))
#     #     # self.assertAlmostEqual(diff, 0)

