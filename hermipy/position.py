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

import numpy as np
import numpy.linalg as la
import sympy as sym


class Position:

    very_small = 1e-10

    @staticmethod
    def tensorize(args):
        dim, mean, cov = 0, [], []
        for a in args:
            assert type(a) is Position
            assert a.is_diag
            dim += a.dim
            mean.extend(a.mean)
            cov.extend(np.diag(a.cov))
        return Position(dim=dim, mean=mean, cov=np.diag(cov))

    @staticmethod
    def inner(p1, p2):
        assert p1.is_diag and p2.is_diag
        diff1 = [d for d in p1.dirs if d not in p2.dirs]
        diff2 = [d for d in p2.dirs if d not in p1.dirs]
        dirs_result = sorted(diff1 + diff2)
        dim, mean, cov = len(dirs_result), [], []
        for d in dirs_result:
            pos = p1 if d in p1.dirs else p2
            index = pos.dirs.index(d)
            mean.append(pos.mean[index])
            cov.append(pos.cov[index][index])
        return Position(dim=dim, mean=mean, cov=np.diag(cov), dirs=dirs_result)

    def __init__(self, dim=None, mean=None, cov=None, dirs=None):

        if mean is not None:
            self.dim = len(mean)
        elif cov is not None:
            self.dim = len(cov)
        else:
            assert dim is not None
            self.dim = dim

        self.mean = np.zeros(self.dim) if mean is None \
            else np.asarray(mean, float)

        self.cov = np.eye(self.dim) if cov is None \
            else np.asarray(cov, float)

        self.dirs = list(range(dim)) if dirs is None else dirs

        eigval, eigvec = la.eig(self.cov)
        self.factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

        diag_cov = np.diag(np.diag(self.cov))
        self.is_diag = la.norm(self.cov - diag_cov, 2) < 1e-10

    def __eq__(self, other):
        return self.dim == other.dim \
            and la.norm(self.mean - other.mean, 2) < self.very_small \
            and la.norm(self.cov - other.cov, 2) < self.very_small

    def __mul__(self, other):
        assert type(other) is Position
        return Position.tensorize([self, other])

    def __hash__(self):
        return hash(frozenset({
            self.dim,
            hash(frozenset(self.mean.flatten())),
            hash(frozenset(self.cov.flatten())),
            hash(frozenset(self.dirs))}))

    def weight(self):
        var = [sym.symbols('v' + str(d), real=True) for d in self.dirs]
        inv_cov = la.inv(self.cov)
        potential = 0.5 * inv_cov.dot(var - self.mean).dot(var - self.mean)
        normalization = 1/(np.sqrt((2*np.pi)**self.dim * la.det(self.cov)))
        return normalization * sym.exp(-potential)

    def weights(self):
        assert self.is_diag
        return [self.project([i]).weight() for i in range(self.dim)]

    def project(self, directions):
        if type(directions) is int:
            directions = [directions]
        assert self.is_diag
        dim = len(directions)
        mean = np.zeros(dim)
        cov = np.zeros((dim, dim))
        for i in range(dim):
            d = directions[i]
            assert d < self.dim
            mean[i] = self.mean[d]
            cov[i][i] = self.cov[d][d]
        return Position(dim=dim, mean=mean, cov=cov, dirs=directions)
