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

#  TODO: Add directions in Position (urbain, 06 Jun 2018)

import hermipy.core as core
# import hermipy.lib as lib
import hermipy.position as pos

from scipy.special import binom
import numpy as np
import numpy.linalg as la
import ipdb


very_small = 1e-10


class Series:

    @staticmethod
    def tensorize(args):
        assert len(args) > 1
        vecs = []
        for a in args:
            assert type(a) is Series
            vecs.append(a.coeffs)
        tens_vec = core.tensorize(vecs)
        tens_pos = pos.Position.tensorize([a.position for a in args])
        return Series(tens_vec, tens_pos)

    def __init__(self, coeffs, position,
                 degree=None, norm=False, index_set="triangle"):
        self.coeffs = coeffs/la.norm(coeffs, 2) if norm else coeffs
        self.position = position
        self.index_set = index_set

        if degree is None:
            dim, npolys = self.position.dim, len(self.coeffs)
            self.degree = core.bissect_degree(dim, npolys, index_set=index_set)
        else:
            self.degree = degree

    def __eq__(self, other):
        assert type(other) is Series
        return self.position == other.position \
            and la.norm(self.coeffs - other.coeffs) < very_small

    def __add__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_coeffs = self.coeffs + other

        elif type(other) is Series:
            assert self.position == other.position
            assert self.index_set == other.index_set
            new_coeffs = self.coeffs + other.coeffs

        else:
            raise TypeError("Invalid type!)")

        return Series(new_coeffs, self.position)

    def __mul__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_coeffs = self.coeffs * other
            return Series(new_coeffs, self.position)

        elif type(other) is Series:
            assert self.index_set == other.index_set
            return Series.tensorize([self, other])

        else:
            raise TypeError("Invalid type: " + str(type(other)))

    def inner(self, other):
        assert type(self) is Series and type(other) is Series
        assert self.index_set == other.index_set
        assert self.degree == other.degree
        d1, d2 = self.position.dirs, other.position.dirs
        result = core.inner(self.coeffs, other.coeffs, d1, d2,
                            index_set=self.index_set)
        inner_pos = pos.Position.inner(self.position, other.position)
        return Series(result, inner_pos, degree=self.degree)

    def project(self, directions):
        if type(directions) is int:
            directions = [directions]
        directions = core.to_numeric(directions)
        p_coeffs = core.project(self.coeffs, self.position.dim, directions,
                                index_set=self.index_set)
        p_pos = self.position.project(directions)
        return Series(p_coeffs, p_pos)

    def subdegree(self, degree):
        assert degree <= self.degree
        n_polys = int(binom(degree + self.position.dim, degree))
        coeffs = self.coeffs[0:n_polys]
        return Series(coeffs, self.position, degree=degree,
                      index_set=self.index_set)


    # def to_cross(self, index_set):
    #     list_cross = 
