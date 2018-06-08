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
import matplotlib.pyplot as plt
import hermipy.quad as hm


def plot_hf(degrees, quad_num, quad_visu, ax):
    dim = len(degrees)
    adim_width = [np.sqrt(2) * np.sqrt(2*d + 1) for d in degrees]
    mean, cov, bounds = quad_num.mean, quad_num.cov, []
    for i in range(dim):
        width = adim_width[i] * np.sqrt(cov[i][i])
        bounds.append([mean[i] - width, mean[i] + width])
    if dim >= 1:
        ax.axvline(x=bounds[0][0])
        ax.axvline(x=bounds[0][1])
    if dim == 2:
        ax.axhline(y=bounds[1][0])
        ax.axhline(y=bounds[1][1])
    deg_max = sum(degrees)
    all_indices = hm.multi_indices(dim, deg_max)
    hf = np.zeros(len(all_indices))
    hf[all_indices.index(tuple(degrees))] = 1
    hf = quad_num.series(hf)
    factor = quad_num.factor_mapping()
    return quad_visu.plot(hf, deg_max, factor, ax)


def plot_projections(series, quad, factors, degree):
    degrees = np.arange(degree + 1)
    series_x = series.project('x')
    series_y = series.project('y')
    quad_x = quad.project('x')
    quad_y = quad.project('y')
    fig, ((ax11, ax12), (ax21, ax22)) = plt.subplots(2, 2)
    quad_x.plot(series_x, degree, factors['x'], ax11)
    ax11.set_title("Marginal distribution in x direction")
    ax12.bar(degrees, series_x.coeffs)
    ax12.set_title("Hermite coefficients of x marginal")
    quad_y.plot(series_y, degree, factors['y'], ax21)
    ax21.set_title("Marginal distribution in y direction")
    ax22.bar(degrees, series_y.coeffs)
    ax22.set_title("Hermite coefficients of y marginal")
    plt.show()
