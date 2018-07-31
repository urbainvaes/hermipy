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

import hermipy.cache as cache
import hermipy.core as core
import hermipy.lib as lib
import hermipy.position as pos
import hermipy.series as hs
import hermipy.settings as rc
import hermipy.function as func
import hermipy.stats as stats
import scipy.sparse as ss
import numpy as np
import numpy.linalg as la
import scipy.sparse.linalg as las
import matplotlib


very_small = 1e-10


class Varf:

    @staticmethod
    def tensorize(args, sparse=None):

        for a in args:
            assert type(a) is Varf

        if len(args) == 1:
            return args[0]

        sparse = rc.settings['sparse'] if sparse is None else sparse
        index_set, degree = args[0].index_set, args[0].degree
        mats = {}
        factor = 1
        for a in args:
            assert type(a) is Varf
            assert a.index_set == index_set
            assert a.degree == degree
            key = frozenset(a.position.dirs)
            mats[key] = a.matrix
            factor *= a.factor.sym

        matrix = core.tensorize(mats, sparse=sparse, index_set=index_set)
        position = pos.Position.tensorize([a.position for a in args])
        factor = func.Function(factor, dirs=position.dirs)
        return Varf(matrix, position, factor=factor, index_set=index_set)

    def __init__(self, matrix, position, factor=1, index_set="triangle"):
        self.is_sparse = isinstance(matrix, ss.csr_matrix)
        self.matrix = matrix
        self.position = position
        self.index_set = index_set

        dim, npolys = self.position.dim, self.matrix.shape[0]
        self.degree = core.iterator_get_degree(dim, npolys,
                                               index_set=index_set)
        assert len(self.multi_indices()) == self.matrix.shape[0]
        self.factor = func.Function(factor, dirs=self.position.dirs)

    def __eq__(self, other):
        assert type(other) is Varf
        norm_func = las.norm if self.is_sparse and other.is_sparse else la.norm
        return self.position == other.position \
            and norm_func(self.matrix - other.matrix) < very_small

    def __add__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_matrix = self.matrix + other

        elif type(other) is Varf:
            assert self.position == other.position
            assert self.index_set == other.index_set
            assert self.factor == other.factor
            new_matrix = self.matrix + other.matrix

        else:
            raise TypeError("Invalid type!)")

        return Varf(new_matrix, self.position,
                    factor=self.factor, index_set=self.index_set)

    def __mul__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_matrix = self.matrix * other
            return Varf(new_matrix, self.position,
                        factor=self.factor, index_set=self.index_set)

        elif type(other) is Varf:
            assert self.index_set == other.index_set
            return Varf.tensorize([self, other])

        else:
            raise TypeError("Invalid type: " + str(type(other)))

    __rmul__ = __mul__
    __radd__ = __add__

    def __sub__(self, other):
        return self + (other*(-1))

    def __call__(self, series):
        assert self.position == series.position
        assert self.index_set == series.index_set
        assert type(series) is hs.Series
        coeffs = self.matrix.dot(series.coeffs)
        return hs.Series(coeffs, self.position,
                         factor=self.factor, index_set=self.index_set)

    def project(self, directions):
        if type(directions) is int:
            directions = [directions]
        rel_dirs = [self.dirs.index(d) for d in directions]
        p_matrix = core.project(self.matrix, self.position.dim, rel_dirs,
                                index_set=self.index_set)
        p_pos = self.position.project(directions)
        factor = self.factor.project(directions)
        return Varf(p_matrix, p_pos,
                    factor=factor, index_set=self.index_set)

    def subdegree(self, degree):
        assert degree <= self.degree
        n_polys = core.iterator_size(self.position.dim, degree)
        kwargs = {'order': 'C'} if not self.is_sparse else {}
        matrix = self.matrix[0:n_polys, 0:n_polys].copy(**kwargs)
        return Varf(matrix, self.position,
                    factor=self.factor, index_set=self.index_set)

    def multi_indices(self):
        return core.iterator_list_indices(self.position.dim, self.degree,
                                          index_set=self.index_set)

    def to_cross(self, degree):
        assert self.index_set == "triangle"
        assert degree + self.position.dim - 1 <= self.degree
        inds_triangle = lib.cross_in_triangle(self.position.dim, degree)
        matrix = self.matrix[np.ix_(inds_triangle, inds_triangle)]
        return Varf(matrix, self.position,
                    factor=self.factor, index_set="cross")

    def solve(self, series):
        assert self.position == series.position
        assert self.index_set == series.index_set
        solve = cache.cache(quiet=True)(las.spsolve if self.is_sparse
                                        else la.solve)
        solution = solve(self.matrix, series.coeffs)
        return hs.Series(solution, position=self.position,
                         factor=self.factor, index_set=self.index_set)

    @stats.log_stats()
    def eigs(self, **kwargs):
        eigs = cache.cache(quiet=True)(las.eigs)
        eig_vals, eig_vecs = eigs(self.matrix, **kwargs)
        result = []
        for v in eig_vecs.T:
            coeffs = np.real(v)
            series = hs.Series(coeffs, self.position,
                               factor=self.factor, index_set=self.index_set)
            result.append(series)
        return eig_vals, result

    def plot(self, ax=None, lines=True):
        show_plt = ax is None

        if show_plt:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1)

        rows, cols, data = [], [], []
        if self.is_sparse:
            s = self.matrix.tocoo()
            rows, cols, data = s.row, s.col, s.data
        else:
            for i in range(self.matrix.shape[0]):
                for j in range(self.matrix.shape[1]):
                    ij = self.matrix[i][j]
                    if abs(ij) > 1e-9:
                        rows.append(i), cols.append(j), data.append(abs(ij))

        pl = ax.scatter(rows, cols, c=data, cmap='ocean_r',
                        norm=matplotlib.colors.LogNorm())
        ax.set_xlim([0, self.matrix.shape[0]])
        ax.set_ylim([0, self.matrix.shape[1]])

        if lines:
            degree = 0
            for i, m in enumerate(self.multi_indices()):
                if sum(m) > degree:
                    degree = degree + 1
                    ax.axvline(x=i, ymin=0, ymax=i)
                    ax.axhline(y=i, xmin=0, xmax=i)

        ax.invert_yaxis()

        if show_plt:
            plt.colorbar(pl, ax=ax)
            plt.show()
            plt.close()
        else:
            return pl
