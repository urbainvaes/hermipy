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

import hermipy as hm
import numpy as np


very_small = 1e-10


class FourierQuad:

    def __init__(self, n_points, bounds=None,
                 position=None, factor=1, dirs=None):

        # Override bounds, dirs
        if position is not None:
            dim = position.dim
            dirs = position.dirs
            bounds = []
            for i in dim:
                mean, width = position.mean[i], position.cov[i][i]
                bounds.append((mean - width/2, mean + width/2))

        elif bounds is not None:
            dim = len(bounds)
            mean, cov = np.zeros(dim), np.zeros(dim, dim)
            for i in dim:
                mean[i] = (bounds[i][0] + bounds[i][1])/2
                cov[i][i] = bounds[i][1] - bounds[i][1]

            if dirs is not None:
                assert len(dirs) == dim
            else:
                dirs = range(dim)

            position = hm.Position(dim=dim, mean=mean, cov=cov, dirs=dirs)

        dim, self.nodes = len(bounds), []
        for i in range(dim):
            left, right = bounds[i][0], bounds[i][1]
            self.nodes.append(np.linspace(left, right, n_points))

        # cov contains bounds!
        mean = (left/2 + right/2 for (left, right) in bounds)
        cov = np.diag([right - left for (left, right) in bounds])

        self.position = position if position is not None else \
            hm.Position(dim=dim, dirs=dirs, mean=mean, cov=cov)

        assert self.position.is_diag

        self.factor = hm.Function(factor, dirs=self.position.dirs)
        self.factor.sym = self.factor.sym.expand()
