import hermite.quad as hm
import hermite.settings as rc
import hermite.series as series
import hermite.position as position

import unittest
import numpy as np
import numpy.linalg as la

import ipdb

settings = {'cache': False, 'cachedir': '/tmp/test_hermite'}
rc.settings.update(settings)


class TestInner(unittest.TestCase):

    def test_inner_positions(self):
        m1, m2 = [1, 2], [3, 4, 5]
        c1, c2 = np.diag([1, 2]), np.diag([3, 4, 5])
        p1 = position.Position(mean=m1, cov=c1, dirs=[0, 1])
        p2 = position.Position(mean=m2, cov=c2, dirs=[1, 3, 4])
        result = position.Position.inner(p1, p2)
        self.assertTrue(result.dirs == [0, 3, 4])
        self.assertAlmostEqual(la.norm(result.mean - [1, 4, 5]), 0)
        self.assertAlmostEqual(la.norm(result.cov - np.diag([1, 4, 5])), 0)

    def test_simple(self):
        n_points, degree = 200, 10
        fx, fy = 'exp(x)', '1'
        quad_xy = hm.Quad.gauss_hermite(n_points, dim=2, dirs=[0, 1])
        quad_y = hm.Quad.gauss_hermite(n_points, dim=1, dirs=[1])
        series_xy = quad_xy.transform(fx, degree)
        series_y = quad_y.transform(fy, degree)
        inner = series.Series.inner(series_xy, series_y)
        self.assertTrue(series_xy.position.dirs == [0, 1])
        self.assertTrue(series_y.position.dirs == [1])
        self.assertTrue(inner.position.dirs == [0])
        projection = series_xy.project(0)
        self.assertTrue(inner == projection)
