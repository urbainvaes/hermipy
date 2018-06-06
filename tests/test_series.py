import hermite.quad as hm
import hermite.settings as rc
import hermite.series as hs
import hermite.position as position

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

import ipdb

settings = {'cache': False, 'cachedir': '/tmp/test_hermite'}
rc.settings.update(settings)
if not os.path.exists(settings['cachedir']):
    os.makedirs(settings['cachedir'])


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

#     def test_simple(self):
#         n_points, degree = 200, 10
#         function = 'exp(x)'
#         quad_1d = hm.Quad.gauss_hermite(n_points, dim=1)
#         quad_2d = hm.Quad.gauss_hermite(n_points, dim=2)
#         series_1 = quad_1d.transform(function, degree)
#         series_2 = quad_2d.transform(function, degree)
#         projection = core.inner()
#         diff = (la.norm(varf_1d - projection, 2))
#         self.assertAlmostEqual(diff, 0)
