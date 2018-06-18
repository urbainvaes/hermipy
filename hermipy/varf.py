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
import hermipy.lib as lib
import hermipy.position as pos

from scipy.special import binom
import scipy.sparse as ss

import numpy as np
import numpy.linalg as la

import ipdb


very_small = 1e-10


class Varf:

    @staticmethod
    def tensorize(args, sparse=False):
        assert len(args) > 1
        mats = []
        for a in args:
            assert type(a) is Varf
            if type(a.matrix) is ss.csr_matrix:
                mats.append(a.matrix.todense().A)
            else:
                mats.append(a.matrix)
        tens_mat = core.tensorize(mats, sparse=sparse)
        tens_pos = pos.Position.tensorize([a.position for a in args])
        return Varf(tens_mat, tens_pos)

    def __init__(self, matrix, position,
                 degree=None, index_set="triangle"):
        self.matrix = matrix
        self.is_sparse = isinstance(matrix, ss.csr_matrix)
        self.position = position
        self.index_set = index_set

        if degree is None:
            dim, npolys = self.position.dim, self.matrix.shape[0]
            self.degree = core.bissect_degree(dim, npolys, index_set=index_set)
        else:
            self.degree = degree

    def __eq__(self, other):
        assert type(other) is Varf
        return self.position == other.position \
            and la.norm(self.matrix - other.matrix) < very_small

    def __add__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_matrix = self.matrix + other

        elif type(other) is Varf:
            assert self.position == other.position
            assert self.index_set == other.index_set
            new_matrix = self.matrix + other.matrix

        else:
            raise TypeError("Invalid type!)")

        return Varf(new_matrix, self.position)

    def __mul__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_matrix = self.matrix * other
            return Varf(new_matrix, self.position)

        elif type(other) is Varf:
            assert self.index_set == other.index_set
            return Varf.tensorize([self, other])

        else:
            raise TypeError("Invalid type: " + str(type(other)))

    def project(self, directions):
        if type(directions) is int:
            directions = [directions]
        directions = core.to_numeric(directions)
        p_matrix = core.project(self.matrix, self.position.dim, directions,
                                index_set=self.index_set)
        p_pos = self.position.project(directions)
        return Varf(p_matrix, p_pos)

    def subdegree(self, degree):
        assert degree <= self.degree
        n_polys = int(binom(degree + self.position.dim, degree))
        matrix = self.matrix[0:n_polys][0:n_polys]
        return Varf(matrix, self.position, degree=degree,
                    index_set=self.index_set)
