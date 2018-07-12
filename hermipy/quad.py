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
import hermipy.settings as rc
import hermipy.function as symfunc
import hermipy.series as hs
import hermipy.varf as hv
import hermipy.position as pos
import hermipy.stats as stats

import numpy as np
import numpy.linalg as la
import sympy as sym


very_small = 1e-10


class Quad:

    def __init__(self, nodes, weights,
                 mean=None, cov=None, dirs=None, position=None):

        """Create a quadrature object

        Args:
            nodes: The nodes of the quadrature.
            weights: The weights of the quadrature, in the basis.

        Returns:
            True if successful, False otherwise.

        Note:
            Do not include the  parameter in the Args section.

        """

        self.nodes = np.asarray(nodes, float)
        self.weights = np.asarray(weights, float)

        dim = len(self.nodes)
        self.position = position if position is not None else \
            pos.Position(dim=dim, mean=mean, cov=cov, dirs=dirs)

    @classmethod
    def gauss_hermite(cls, n_points, dim=None, dirs=None, **kwargs):

        if dirs is not None:
            dim = len(dirs)

        if dim is not None:
            n_points = np.full(dim, n_points)

        elif isinstance(n_points, int):
            n_points = [n_points]

        nodes, weights = lib.hermegauss_nd(n_points)
        return cls(nodes, weights, dirs=dirs, **kwargs)

    @classmethod
    def newton_cotes(cls, n_points, extrema, **kwargs):
        nodes, weights = [], []
        for i in range(len(extrema)):
            nodes.append(np.linspace(-extrema[i], extrema[i], n_points[i]))
            mesh_size = 2*extrema[i]/(n_points[i] - 1)
            weights_simpson = np.zeros(n_points[i]) + 1
            weights_simpson[0], weights_simpson[-1] = .5, .5
            gaussian_weight = 1/np.sqrt(2*np.pi) * np.exp(-nodes[-1]**2/2.)
            weights.append(weights_simpson * gaussian_weight * mesh_size)
        return cls(nodes, weights, **kwargs)

    def __mul__(self, other):
        assert self.position.is_diag and other.position.is_diag
        nodes = [*self.nodes, *other.nodes]
        weights = [*self.weights, *other.weights]
        position = self.position * other.position
        return Quad(nodes, weights, position=position)

    def __eq__(self, other):
        assert type(other) is Quad

        return self.position == other.position \
            and la.norm(self.nodes - other.nodes) < very_small \
            and la.norm(self.weights - other.weights) < very_small

    def hash_quad(argument):
        if isinstance(argument, Quad):
            return hash(argument)
        raise ValueError("Argument type not supported")

    def __hash__(self):
        return self.hash

    def tensorize_at(arg_num):
        def tensorize_arg(func):
            def wrapper(*args, **kwargs):

                do_tensorize = rc.settings['tensorize']
                if 'tensorize' in kwargs:
                    do_tensorize = kwargs['tensorize']
                    del kwargs['tensorize']
                if not do_tensorize:
                    return func(*args, **kwargs)

                quad, function = args[0], args[arg_num]
                if not quad.position.is_diag or \
                        isinstance(function, np.ndarray) or \
                        quad.position.dim is 1:
                    return func(*args, **kwargs)

                results = []
                if not isinstance(function, symfunc.Function):
                    function = symfunc.Function(function,
                                                dirs=quad.position.dirs)

                for add in function.split():
                    if len(add) == 2:
                        new_args = list(args).copy()
                        new_args[arg_num] = add[0]
                        results.append(func(*new_args, **kwargs)*float(add[1]))
                        continue

                    func_dirs = []
                    for d in range(quad.position.dim):
                        new_args = list(args).copy()
                        new_args[0] = quad.project(d)
                        new_args[arg_num] = add[d]
                        func_dir = func(*new_args, **kwargs)
                        func_dirs.append(func_dir)
                    if rc.settings['debug']:
                        print("Tensorizing results")
                    kwargs_func = {'sparse': kwargs['sparse']} \
                        if 'sparse' in kwargs else {}
                    t = type(func_dirs[0])
                    tens_fun = t.tensorize if t is hv.Varf or t is hs.Series \
                        else core.tensorize
                    tensorized = tens_fun(func_dirs, **kwargs_func)
                    results.append(tensorized*float(add[-1]))

                return sum(results[1:], results[0])
            return wrapper
        return tensorize_arg

    def mapped_nodes(self):
        coords_nodes = []
        for i in range(self.position.dim):
            coord = 'v[{}]'.format(i)
            coords_nodes.append(self.discretize(coord))
        return np.asarray(np.vstack(coords_nodes)).T

    def discretize(self, f):
        if not isinstance(f, symfunc.Function):
            f = symfunc.Function(f, dirs=self.position.dirs)
        f = f.as_string(format='array', toC=True)
        function = core.discretize(f, self.nodes,
                                   self.position.mean, self.position.factor)
        return function

    @tensorize_at(1)
    def integrate(self, f_grid, l2=False):
        if not isinstance(f_grid, np.ndarray):
            f_grid = self.discretize(f_grid)
        if l2:
            w_grid = self.discretize(self.position.weight())
            f_grid = f_grid / w_grid
        return core.integrate(f_grid, self.nodes, self.weights)

    # Norm 1 or 2, in weighted or not
    def norm(self, function, n=2, l2=False):
        if n is 2:
            return np.sqrt(self.integrate(function**2, l2=l2))
        elif n is 1:
            return self.integrate(abs(function), l2=l2)

    def transform(self, f_grid, degree, norm=False,
                  index_set="triangle", significant=0):
        if not isinstance(f_grid, np.ndarray):
            f_grid = self.discretize(f_grid)
        coeffs = core.transform(degree, f_grid, self.nodes, self.weights,
                                forward=True, index_set=index_set)
        return hs.Series(coeffs, self.position, norm=norm,
                         index_set=index_set, significant=significant)

    def eval(self, series):
        if type(series) is np.ndarray:
            series = hs.Series(series, self.position)
        #  FIXME: Only orientation, not positions
        # assert self.position == series.position
        degree, coeffs = series.degree, series.coeffs
        inv = la.inv(series.position.factor)
        translation = inv.dot(self.position.mean - series.position.mean)
        factor = inv * self.position.factor
        if la.norm(factor - np.diag(np.diag(factor)), 2) > 1e-8:
            raise ValueError("Incompatible covariance matrices")
        mapped_nodes = self.nodes.copy()
        for i in range(len(self.nodes)):
            mapped_nodes[i] = self.nodes[i] * factor[i][i] + translation[i]
        return core.transform(degree, coeffs, mapped_nodes,
                              self.weights, forward=False,
                              index_set=series.index_set)

    @tensorize_at(1)
    @stats.debug()
    @stats.log_stats()
    def varf(self, f_grid, degree, sparse=False, index_set="triangle"):
        if not isinstance(f_grid, np.ndarray):
            f_grid = self.discretize(f_grid)
        var = core.varf(degree, f_grid, self.nodes, self.weights,
                        sparse=sparse, index_set=index_set)
        return hv.Varf(var, self.position, index_set=index_set)

    @stats.debug()
    @stats.log_stats()
    def varfd(self, function, degree, directions, sparse=False,
              index_set="triangle"):
        var = self.varf(function, degree, sparse=sparse, index_set=index_set)
        mat = var.matrix
        eigval, _ = la.eig(self.position.cov)
        for d in directions:
            mat = core.varfd(self.position.dim, degree, d, mat,
                             sparse=sparse, index_set=index_set)
            mat = mat/np.sqrt(eigval[d])
        return hv.Varf(mat, self.position, index_set=index_set)

    @stats.debug()
    @stats.log_stats()
    def discretize_op(self, op, func, degree, order,
                      sparse=None, index_set="triangle"):

        assert len(func.args) <= self.position.dim
        sparse = rc.settings['sparse'] if sparse is None else sparse
        mat_operator = 0.
        splitop, mult = lib.split_operator(op, func, order)
        for m, coeff in zip(mult, splitop):
            d_vector = sum([[i]*m[i] for i in range(self.position.dim)], [])
            varf_part = self.varfd(coeff, degree, d_vector, sparse=sparse,
                                   index_set=index_set)
            mat_operator = varf_part + mat_operator
        return mat_operator
    # Only works with ints
    def project(self, directions):
        """Project the quadrature.

        Args:
            directions: a list of integers containing the directions along
            which the quadrature is to be projected.

        Returns:
            A `Quad` object.

        """
        if type(directions) is int:
            directions = [directions]
        dim = len(directions)
        nodes, weights = [], []
        for i in range(dim):
            d = directions[i]
            assert d < self.position.dim
            nodes.append(self.nodes[d])
            weights.append(self.weights[d])
        pos = self.position.project(directions)
        return Quad(nodes, weights, position=pos)

    def series(self, coeffs, norm=False, index_set="triangle"):
        """Return a series from a vector of coefficients.

        Args:
            coeffs: a numpy array of coefficients.

            norm: whether or not to normalize the coefficients (using the
                  little L² norm).

            index_set: the index set corresponding to the series.

        Returns:
            A `Series` object.

        """
        return hs.Series(coeffs, self.position, norm=norm,
                         index_set=index_set)

    def plot(self, arg, factor=None, ax=None, bounds=True):
        assert self.position.is_diag

        show_plt = ax is None
        if show_plt:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1)

        if isinstance(arg, tuple(sym.core.all_classes)):
            arg = symfunc.Function(arg, dirs=self.position.dirs)

        if type(arg) is symfunc.Function:
            assert factor is None
            solution = self.discretize(arg)

        elif type(arg) is hs.Series:

            if factor is None:
                factor = self.position.weight()

            if not isinstance(factor, np.ndarray):
                factor = symfunc.Function(factor, dirs=self.position.dirs)
                factor = self.discretize(factor)

            series = arg
            solution = self.eval(series)*factor

            if bounds:

                bounds, adim_width = [], np.sqrt(2)*np.sqrt(2*series.degree+1)
                for i in range(series.position.dim):
                    mean = series.position.mean[i]
                    cov = series.position.cov[i][i]
                    bounds.append([mean - adim_width * np.sqrt(cov),
                                   mean + adim_width * np.sqrt(cov)])

                if self.position.dim >= 1:
                    ax.axvline(x=bounds[0][0])
                    ax.axvline(x=bounds[0][1])
                if self.position.dim == 2:
                    ax.axhline(y=bounds[1][0])
                    ax.axhline(y=bounds[1][1])

        else:
            raise TypeError("Unsupported type: " + str(type(arg)))

        n_nodes = []
        r_nodes = []
        for i in range(self.position.dim):
            direction = symfunc.Function.xyz[self.position.dirs[i]]
            n_nodes.append(len(self.nodes[i]))
            r_nodes.append(self.project(i).discretize(direction))
        solution = solution.reshape(*n_nodes).T

        if self.position.dim == 1:
            plot = ax.plot(*r_nodes, solution)
        elif self.position.dim == 2:
            plot = ax.contourf(*r_nodes, solution, 100)

        min, max = np.min(solution), np.max(solution)
        ax.set_title("Min: {:.3f}, Max: {:.3f}".format(min, max))

        if show_plt:
            if self.position.dim == 2:
                plt.colorbar(plot, ax=ax)
            plt.show()
        else:
            return plot

    def streamlines(self, fx, fy, factor=None, ax=None, **kwargs):
        assert self.position.is_diag
        assert type(fx) is type(fy)
        assert self.position.dim is 2

        show_plt = ax is None
        if show_plt:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1)

        if isinstance(fx, tuple(sym.core.all_classes)):
            fx = symfunc.Function(fx, dirs=self.position.dirs)
            fy = symfunc.Function(fy, dirs=self.position.dirs)

        if type(fx) is symfunc.Function:
            assert factor is None
            fx, fy = self.discretize(fx), self.discretize(fy)

        elif type(fx) is hs.Series:

            if factor is None:
                factor = self.position.weight()

            if not isinstance(factor, np.ndarray):
                factor = symfunc.Function(factor, dirs=self.position.dirs)
                factor = self.discretize(factor)

            sx, sy = self.eval(fx)*factor, self.eval(fy)*factor

        else:
            raise TypeError("Unsupported type: " + str(type(fx)))

        n_nodes = []
        r_nodes = []
        for i in range(self.position.dim):
            direction = symfunc.Function.xyz[self.position.dirs[i]]
            n_nodes.append(len(self.nodes[i]))
            r_nodes.append(self.project(i).discretize(direction))
        sx, sy = sx.reshape(*n_nodes).T, sy.reshape(*n_nodes).T
        # magnitude = sx**2 + sy**2
        x, y = r_nodes[0], r_nodes[1]

        # See https://github.com/matplotlib/matplotlib/issues/9269
        def calc_psi2d(x, y, fx, fy):  # solenoidal flows
            # psi by integrating dpsi = -fy*dx + fx*dy, psi[0,0]=0
            ny, nx = fx.shape
            psi = np.zeros((ny, nx))
            for jx in range(1, nx):
                psi[0, jx] = psi[0, jx-1] - fy[0, jx] * (x[jx] - x[jx-1])
            for jy in range(1, ny):
                psi[jy, :] = psi[jy-1, :] + fx[jy, :] * (y[jy] - y[jy-1])
            return psi

        psi = calc_psi2d(x, y, sx, sy)
        min, max, thresh = np.min(psi), np.max(psi), .1
        if abs(min) < thresh*(max - min):
            min = thresh*(max - min)
        if abs(max) < thresh*(max - min):
            max = - thresh*(max - min)
        levels = np.linspace(min, max, 10)
        streams = ax.contour(x, y, psi, levels=levels, **kwargs)

        # ax.quiver(x, y, sx, sy)
        # streams = ax.streamplot(r_nodes[0], r_nodes[1], sx, sy,
        #                         color=magnitude, density=0.6, cmap='autumn')

        return plt.show() if show_plt else streams
