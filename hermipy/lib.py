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

import numpy.polynomial.hermite_e as herm
import numpy as np
import hermipy.core as core


def hermegauss_nd(n_points):
    dim = len(n_points)
    nodes_multidim = []
    weights_multidim = []
    for i in range(dim):
        nodes_1d, weights_1d = herm.hermegauss(n_points[i])
        weights_1d = weights_1d/np.sqrt(2*np.pi)  # Normalize
        nodes_multidim.append(nodes_1d)
        weights_multidim.append(weights_1d)
    return nodes_multidim, weights_multidim


def cross_in_triangle(dim, degree):
        list_cross = core.iterator_list_indices(dim, degree, index_set="cross")
        return [core.iterator_index(m) for m in list_cross]


def finest_common(set1, set2):

    set1, set2 = set1.copy(), set2.copy()

    difference = set1 - set2
    if difference == set():
        return set1

    s = list(difference)[0]
    set1.remove(s)

    while (any(s1.intersection(s) != set() for s1 in set1) or
           any(s2.intersection(s) != set() for s2 in set2)):

        for s2 in list(set2):
            if s.intersection(s2) != set():
                s = s.union(s2)
                set2.remove(s2)

        for s1 in list(set1):
            if s.intersection(s1) != set():
                s = s.union(s1)
                set1.remove(s1)

    set1.add(s)
    set2.add(s)

    return finest_common(set1, set2)
