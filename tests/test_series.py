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

import hermipy.quad as hm
import hermipy.settings as rc
import hermipy.series as series
import hermipy.position as position
import unittest
import numpy as np
import numpy.linalg as la


class TestSeries(unittest.TestCase):

    def setUp(self):

        settings = {
            'cache': False,
            'cachedir': '/tmp/test_hermite',
            'tensorize': False,
            'sparse': False,
            'debug': False,
            }

        rc.settings.update(settings)

    def test_inner_positions(self):
        m1, m2 = [1, 2], [3, 4, 5]
        c1, c2 = np.diag([1, 2]), np.diag([3, 4, 5])
        p1 = position.Position(mean=m1, cov=c1, dirs=[0, 1])
        p2 = position.Position(mean=m2, cov=c2, dirs=[1, 3, 4])
        result = position.Position.inner(p1, p2)
        self.assertTrue(result.dirs == [0, 3, 4])
        self.assertAlmostEqual(la.norm(result.mean - [1, 4, 5]), 0)
        self.assertAlmostEqual(la.norm(result.cov - np.diag([1, 4, 5])), 0)

    def aux_test_simple(self, index_set):
        # import ipdb; ipdb.set_trace()
        n_points, degree = 20, 10
        fx, fy = 'exp(x)', '1'
        quad_xy = hm.Quad.gauss_hermite(n_points, dim=2, dirs=[0, 1])
        quad_y = hm.Quad.gauss_hermite(n_points, dim=1, dirs=[1])
        series_xy = quad_xy.transform(fx, degree, index_set=index_set)
        series_y = quad_y.transform(fy, degree, index_set=index_set)
        inner = series.Series.inner(series_xy, series_y)
        self.assertTrue(series_xy.position.dirs == [0, 1])
        self.assertTrue(series_y.position.dirs == [1])
        self.assertTrue(inner.position.dirs == [0])
        projection = series_xy.project(0)
        self.assertTrue(inner == projection)

    def test_simple_triangle(self):
        self.aux_test_simple("triangle")

    def test_simple_cross(self):
        self.aux_test_simple("cross")

    def test_simple_cube(self):
        self.aux_test_simple("cube")

    def test_consistence_series_1d(self):
        n_points, degree = 200, 40
        fx = 'exp(x - .5)'
        quad = hm.Quad.gauss_hermite(n_points, dim=1)
        series_cross = quad.transform(fx, degree, index_set="cross")
        series_triangle = quad.transform(fx, degree, index_set="triangle")
        diff = la.norm(series_cross.coeffs - series_triangle.coeffs)
        self.assertAlmostEqual(diff, 0.)

    def test_to_cross(self):
        n_points, degree = 200, 50
        fxy = 'exp(x) * cos(y) * y**4'
        quad = hm.Quad.gauss_hermite(n_points, dim=2, dirs=[0, 1])
        series_cross = quad.transform(fxy, degree - 1, index_set="cross")
        series_triangle = quad.transform(fxy, degree, index_set="triangle")
        triangle_to_cross = series_triangle.to_cross(degree - 1)
        self.assertTrue(series_cross == triangle_to_cross)

    def test_approximation_2d(self):
        n_points, degree = 200, 100
        fxy = 'exp(x) * cos(y) * y**4'
        quad = hm.Quad.gauss_hermite(n_points, dim=2, dirs=[0, 1])
        series_cross = quad.transform(fxy, degree, index_set="cross")
        series_triangle = quad.transform(fxy, degree + 1, index_set="triangle")
        exact = quad.discretize(fxy)
        eval_cross = quad.eval(series_cross)
        eval_triangle = quad.eval(series_triangle)
        error_cross = quad.norm(exact - eval_cross)
        error_triangle = quad.norm(exact - eval_triangle)
        self.assertTrue(error_triangle < 1e-10)
        self.assertTrue(error_cross < 1e-3)
