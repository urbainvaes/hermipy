# Todos {{{
# TODO: Implement composite quadrature
# TODO: Ensure directions match
# TODO: Can varfd be tensorized?
# TODO: Improve separability: for example f(x,y) * g(z)
# TODO: Implement composite quadrature (urbain, 02 May 2018)
# TODO: Improve linearize to notice constant * operator
# }}}
import hermite.core as core
import hermite.lib as lib
import hermite.settings as rc
import hermite.function as symfunc
import hermite.series as hs
import hermite.varf as hv

from hermite.cache import cache
from scipy.special import binom

import numpy as np
import numpy.linalg as la
import sympy as sym
import ipdb


very_small = 1e-10


class Quad:

    def __init__(self, nodes, weights, mean=None, cov=None):
        self.nodes = np.asarray(nodes, float)
        self.weights = np.asarray(weights, float)

        self.dim = len(self.nodes)
        self.mean = np.zeros(self.dim) if mean is None \
            else np.asarray(mean, float)
        self.cov = np.eye(self.dim) if cov is None \
            else np.asarray(cov, float)

        eigval, eigvec = la.eig(self.cov)
        self.factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

        diag_cov = np.diag(np.diag(self.cov))
        self.is_diag = la.norm(self.cov - diag_cov, 2) < 1e-10

        self.hash = hash(frozenset({
            self.dim,
            hash(frozenset(self.nodes.flatten())),
            hash(frozenset(self.weights.flatten())),
            hash(frozenset(self.mean.flatten())),
            hash(frozenset(self.cov.flatten()))}))

    def __mul__(self, other):
        assert self.is_diag and other.is_diag
        dim = self.dim + other.dim
        nodes = [*self.nodes, *other.nodes]
        weights = [*self.weights, *other.weights]
        mean = np.zeros(dim)
        cov = np.zeros((dim, dim))
        for i in range(self.dim):
            mean[i] = self.mean[i]
            cov[i][i] = self.cov[i][i]
        for i in range(other.dim):
            off = self.dim
            mean[off + i] = other.mean[i]
            cov[off + i][off + i] = other.cov[i][i]
        return Quad(nodes, weights, mean, cov)

    def __eq__(self, other):
        assert type(other) is Quad

        return self.dim == other.dim \
            and la.norm(self.mean - other.mean, 2) < very_small \
            and la.norm(self.cov - other.cov, 2) < very_small \
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
                if not quad.is_diag or isinstance(function, np.ndarray):
                    return func(*args, **kwargs)

                results = []
                # pdb.set_trace()
                if not isinstance(function, symfunc.Function):
                    function = symfunc.Function(function, dim=quad.dim)

                for add in function.split():
                    if len(add) == 2:
                        new_args = list(args).copy()
                        new_args[arg_num] = add[0]
                        results.append(func(*new_args, **kwargs)*float(add[1]))
                        continue

                    func_dirs = []
                    for d in range(quad.dim):
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
                    # pdb.set_trace()
                    results.append(tensorized*float(add[-1]))

                return sum(results[1:], results[0])
            return wrapper
        return tensorize_arg

    @classmethod
    def gauss_hermite(cls, n_points, dim=None, mean=None, cov=None):
        if dim is not None:
            n_points = np.full(dim, n_points)
        elif isinstance(n_points, int):
            n_points = [n_points]
        nodes, weights = lib.hermegauss_nd(n_points)
        return cls(nodes, weights, mean=mean, cov=cov)

    @classmethod
    def newton_cotes(cls, n_points, extrema, mean=None, cov=None):
        nodes, weights = [], []
        for i in range(len(extrema)):
            nodes.append(np.linspace(-extrema[i], extrema[i], n_points[i]))
            mesh_size = 2*extrema[i]/(n_points[i] - 1)
            weights_simpson = np.zeros(n_points[i]) + 1
            weights_simpson[0], weights_simpson[-1] = .5, .5
            gaussian_weight = 1/np.sqrt(2*np.pi) * np.exp(-nodes[-1]**2/2.)
            weights.append(weights_simpson * gaussian_weight * mesh_size)
        return cls(nodes, weights, mean=mean, cov=cov)

    def mapped_nodes(self):
        coords_nodes = []
        for i in range(self.dim):
            coord = 'v[{}]'.format(i)
            coords_nodes.append(self.discretize(coord))
        return np.asarray(np.vstack(coords_nodes)).T

    def discretize(self, f):
        if not isinstance(f, symfunc.Function):
            f = symfunc.Function(f)
        function = core.discretize(str(f), self.nodes, self.mean, self.factor)
        return function

    @tensorize_at(1)
    def integrate(self, f_grid, l2=False):
        if not isinstance(f_grid, np.ndarray):
            f_grid = self.discretize(f_grid)
        if l2:
            w_grid = self.discretize(self.weight())
            f_grid = f_grid / w_grid
        return core.integrate(f_grid, self.nodes, self.weights)

    # Norm 1 or 2, in weighted or not
    def norm(self, function, n=2, l2=False):
        if n is 2:
            return np.sqrt(self.integrate(function**2, l2=l2))
        elif n is 1:
            return self.integrate(abs(function), l2=l2)

    def transform(self, f_grid, degree, norm=False):
        if not isinstance(f_grid, np.ndarray):
            f_grid = self.discretize(f_grid)
        coeffs = core.transform(degree, f_grid, self.nodes,
                                self.weights, forward=True)
        return hs.Series(coeffs, self.dim, self.mean, self.cov,
                         norm=norm, degree=degree)

    def eval(self, series):
        if type(series) is np.ndarray:
            series = hs.Series(series, self.dim, self.mean, self.cov)
        degree, coeffs = series.degree, series.coeffs
        inv = la.inv(series.factor)
        translation = inv.dot(self.mean - series.mean)
        factor = inv * self.factor
        if la.norm(factor - np.diag(np.diag(factor)), 2) > 1e-8:
            raise ValueError("Incompatible covariance matrices")
        mapped_nodes = self.nodes.copy()
        for i in range(len(self.nodes)):
            mapped_nodes[i] = self.nodes[i] * factor[i][i] + translation[i]
        return core.transform(degree, coeffs, mapped_nodes,
                              self.weights, forward=False)

    @tensorize_at(1)
    def varf(self, f_grid, degree, sparse=False, numpy=True):
        if rc.settings['debug']:
            print("Entering body of Quad.varf")
        if not isinstance(f_grid, np.ndarray):
            f_grid = self.discretize(f_grid)
        var = core.varf(degree, f_grid, self.nodes,
                        self.weights, sparse=sparse)
        return hv.Varf(var, self.dim, self.mean, self.cov, degree=degree)

    def varfd(self, function, degree, directions, sparse=False):
        directions = core.to_numeric(directions)
        var = self.varf(function, degree, sparse=sparse)
        mat = var.matrix
        eigval, _ = la.eig(self.cov)
        for d in directions:
            # ipdb.set_trace()
            mat = core.varfd(self.dim, degree, d, mat)
            mat = mat/np.sqrt(eigval[d])
        return hv.Varf(mat, self.dim, self.mean, self.cov, degree=degree)

    @cache(hash_extend=hash_quad)
    def discretize_op(self, op, func, degree, order, sparse=False):
        npolys = int(binom(degree + self.dim, degree))
        mat_operator = np.zeros((npolys, npolys))
        mult = list(core.multi_indices(self.dim, order))
        splitop = lib.split_operator(op, func, order)
        v = ['x', 'y', 'z']
        for m, coeff in zip(mult, splitop):
            d_vector = sum([[v[i]]*m[i] for i in range(self.dim)], [])
            # ipdb.set_trace()
            mat_operator += self.varfd(coeff, degree, d_vector, sparse=sparse)
        return mat_operator

    #  TODO: Ensure order is right (urbain, Tue 01 May 2018)
    def plot(self, series, factor, ax=None):
        assert self.is_diag
        if not isinstance(factor, np.ndarray):
            factor = self.discretize(factor)
        n_nodes = []
        r_nodes = []
        for i in range(self.dim):
            n_nodes.append(len(self.nodes[i]))
            r_nodes.append(self.project(i).discretize('x'))
        solution = self.eval(series)*factor
        solution = solution.reshape(*n_nodes).T
        if self.dim == 1:
            return ax.plot(*r_nodes, solution)
        elif self.dim == 2:
            return ax.contourf(*r_nodes, solution, 100)

    def weight(self):
        var = [sym.symbols('v' + str(i), real=True) for i in range(self.dim)]
        inv_cov = la.inv(self.cov)
        potential = 0.5 * inv_cov.dot(var - self.mean).dot(var - self.mean)
        normalization = 1/(np.sqrt((2*np.pi)**self.dim * la.det(self.cov)))
        return normalization * sym.exp(-potential)

    #  Only works with ints
    def project(self, directions):
        assert self.is_diag
        if type(directions) is int:
            directions = [directions]
        dim = len(directions)
        nodes, weights = [], []
        mean = np.zeros(dim)
        cov = np.zeros((dim, dim))
        for i in range(len(directions)):
            d = directions[i]
            assert d < self.dim
            nodes.append(self.nodes[d])
            weights.append(self.weights[d])
            mean[i] = self.mean[d]
            cov[i][i] = self.cov[d][d]
        return Quad(nodes, weights, mean, cov)

    def series(self, coeffs, degree=None, norm=False):
        return hs.Series(coeffs,
                         dim=self.dim,
                         mean=self.mean,
                         cov=self.cov,
                         degree=degree,
                         norm=norm)
