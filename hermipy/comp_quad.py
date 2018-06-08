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

#  Composite quadrature {{{


class CompQuad:
    def __init__(self, quads, weights):
        self.quads = quads
        self.weights = weights

    # def integrate(f):

    # def eval(self, degree, nodes):
    #     return eval_simple_quad(self.coeffs, degree, nodes)

# def herm_to_poly(c):
#     herme_coeffs = c/np.sqrt(np.sqrt(2*np.pi)*np.arange(len(c)))
#     return herm.herme2poly(herme_coeffs)
# }}}
# vim: foldmethod=marker foldnestmax=2
