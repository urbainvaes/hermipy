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
import hermipy.function as func

import numpy as np
import numpy.linalg as la

from matplotlib.ticker import MaxNLocator


very_small = 1e-10


class Series:

    @staticmethod
    def tensorize(args):
        assert len(args) > 0 and type(args[0]) is Series
        index_set, degree = args[0].index_set, args[0].degree
        vecs = {}
        for a in args:
            assert type(a) is Series
            assert a.index_set == index_set and a.degree == degree
            key = frozenset(a.position.dirs)
            vecs[key] = a.coeffs
        tens_vec = core.tensorize(vecs, index_set=index_set)
        tens_pos = pos.Position.tensorize([a.position for a in args])
        return Series(tens_vec, tens_pos, index_set=index_set)

    def __init__(self, coeffs, position, norm=False,
                 index_set="triangle", significant=0):
        self.coeffs = coeffs/la.norm(coeffs, 2) if norm else coeffs
        self.position = position
        self.index_set = index_set

        dim, npolys = self.position.dim, len(self.coeffs)
        self.degree = core.iterator_get_degree(dim, npolys,
                                               index_set=index_set)
        assert len(self.multi_indices()) == len(self.coeffs)

        if significant is not 0:
            for i, c in enumerate(self.coeffs):
                self.coeffs[i] = round(self.coeffs[i], significant)

    def __eq__(self, other):
        assert type(other) is Series
        return self.position == other.position \
            and la.norm(self.coeffs - other.coeffs) < very_small

    def __add__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_coeffs = self.coeffs + other

        elif type(other) is Series:
            assert self.position == other.position

            if self.index_set == other.index_set:
                new_coeffs = self.coeffs + other.coeffs

            #  TODO: Add support for addition of different degrees / index sets

        else:
            raise TypeError("Invalid type!)")

        return Series(new_coeffs, self.position)

    def __mul__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_coeffs = self.coeffs * other
            return Series(new_coeffs, self.position, index_set=self.index_set)

        elif type(other) is Series:
            assert self.index_set == other.index_set
            return Series.tensorize([self, other])

        else:
            raise TypeError("Invalid type: " + str(type(other)))

    __rmul__ = __mul__

    def __sub__(self, other):
        return self + other * (-1)

    def __repr__(self):
        m_list = self.multi_indices()
        assert len(m_list) == len(self.coeffs)
        result = ""
        for m, c in zip(m_list, self.coeffs):
            result += str(m) + ": " + str(c) + "\n"
        return result

    def inner(self, other):
        assert type(self) is Series and type(other) is Series
        assert self.index_set == other.index_set
        assert self.degree == other.degree
        d1, d2 = self.position.dirs, other.position.dirs
        result = core.inner(self.coeffs, other.coeffs, d1, d2,
                            index_set=self.index_set)
        inner_pos = pos.Position.inner(self.position, other.position)
        return Series(result, inner_pos, index_set=self.index_set)

    def project(self, directions):
        if type(directions) is not list:
            directions = [directions]
        rel_dirs = [self.position.dirs.index(d) for d in directions]
        p_coeffs = core.project(self.coeffs, self.position.dim, rel_dirs,
                                index_set=self.index_set)
        p_pos = self.position.project(directions)
        return Series(p_coeffs, p_pos, index_set=self.index_set)

    def subdegree(self, degree):
        assert degree <= self.degree
        n_polys = core.iterator_size(self.position.dim, degree)
        coeffs = self.coeffs[0:n_polys]
        return Series(coeffs, self.position, index_set=self.index_set)

    def multi_indices(self):
        return core.iterator_list_indices(self.position.dim, self.degree,
                                          index_set=self.index_set)

    def plot(self, ax=None):

        show_plt = ax is None
        if show_plt:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1)

        m = self.multi_indices()
        if self.position.dim is 1:
            mx = m[:, 0]
            pl = ax.bar(mx, self.coeffs)
        elif self.position.dim is 2:
            coeffs = self.coeffs
            # coeffs = self.coeffs / max(abs(self.coeffs))
            # coeffs = (abs(coeffs) > 1e-10) * coeffs
            mx, my = m[:, 0], m[:, 1]
            # zoom = 1e3/max(mx)
            # positives = zoom * (coeffs > 0) * coeffs
            # negatives = zoom * (coeffs < 0) * coeffs * (-1)
            zeros = abs(coeffs) < 1e-12
            # ax.scatter(mx, my, s=positives, c='g', marker='o')
            pl = ax.scatter(mx, my, c=abs(coeffs),
                            cmap='ocean_r', s=1e2, edgecolor='k')
            ax.scatter(mx, my, s=zeros, c='r', marker='o')
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))

            # for i, txt in enumerate(coeffs):
            #     ax.annotate("{:.2f}".format(txt), (mx[i] + .1, my[i] + .1), size=8)

        ax.set_title("Coefficients of the Hermite expansion")

        if show_plt:
            if self.position.dim is 2:
                plt.colorbar(pl, ax=ax)
            plt.show()
        else:
            return pl

    def to_function(self):
        assert self.position.is_diag

        def rec_a(i):
            return float(1/np.sqrt(i+1))

        def rec_b(i):
            return - float(np.sqrt(i)/np.sqrt(i+1))

        dirs = self.position.dirs
        hermite_dirs = []

        for i, d in enumerate(dirs):
            var = func.Function.xyz[d]
            h0 = func.Function('1', dirs=dirs)
            μ, σ = self.position.mean[i], self.position.cov[i][i]
            f1 = np.sqrt(σ) * (var - μ)
            h1 = func.Function(f1, dirs=dirs)
            hermite_dirs.append([h0, h1])
            for j in range(1, self.degree):
                hermite_dirs[i].append(hermite_dirs[i][-1]
                                       * hermite_dirs[i][1] * rec_a(j)
                                       + hermite_dirs[i][-2] * rec_b(j))

        result = func.Function('0', dirs=dirs)
        for im, m in enumerate(self.multi_indices()):
            result_m = func.Function(self.coeffs[im], dirs=dirs)
            for i in range(0, self.position.dim):
                result_m *= hermite_dirs[i][m[i]]
            result += result_m
        return result

    def to_cross(self, degree):
        assert self.index_set == "triangle"
        assert degree + self.position.dim - 1 <= self.degree
        coeffs = self.coeffs[lib.cross_in_triangle(self.position.dim, degree)]
        return Series(coeffs, self.position, index_set="cross")
