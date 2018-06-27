#!/usr/bin/env python
#
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

import hermipy.equations as eq
import hermipy.quad as quad
import hermipy.series as series

import numpy as np
import numpy.linalg as la

backward = eq.Langevin.backward(1)
forward = eq.Langevin.forward(1)

degree, nquad = 10, 20

σx, σy = 1, 1
quad_num = quad.Quad.gauss_hermite(nquad, dim=2, mean=[0, 0],
                                   cov=[[σx, 0], [0, σy]])

op = quad_num.discretize_op(backward, eq.Langevin.f, degree, 2,
                            sparse=False, index_set="triangle")

rhs = quad_num.transform('y', degree)

# The matrix is singular, so we remove first row and column
solution = la.solve(op.matrix[1:, 1:], rhs.coeffs[1:])
solution = np.array([0, *solution])
sol_series = series.Series(solution, rhs.position, significant=10)

symbolic = sol_series.to_sympy()
