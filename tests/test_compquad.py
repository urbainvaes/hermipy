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

import hermipy.comp_quad as hm
import hermipy.settings as rc
import hermipy.lib as lib

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


class TestCompQuad(unittest.TestCase):

    def setUp(self):

        settings = {
            'cache': False,
            'cachedir': '/tmp/test_hermite',
            'tensorize': False,
            'sparse': False,
            'trails': True,
            'debug': False,
            }

        rc.settings.update(settings)

