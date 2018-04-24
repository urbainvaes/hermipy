#  TODO: Ensure directions match (urbain, Wed 28 Mar 2018 11:55:08 AM BST)
#  TODO: Add degree to series attributes


from .cpp import hermite_cpp as hm
import numpy as np
import inspect
import itertools

import numpy.linalg as la
import numpy.polynomial.hermite_e as herm
import sympy as sym
import re
import math
import matplotlib.pyplot as plt

from scipy.special import binom


def convert(conv_func, *names):
    def convert(function):
        sig = inspect.signature(function)

        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            args_od = ba.arguments
            for key in args_od:
                if key in names:
                    arg = args_od[key]
                    args_od[key] = conv_func(arg)
            return function(**args_od)
        return wrapper
    return convert


def convert_to_cpp_vec(vec):
    cpp_vec = hm.double_vec()
    cpp_vec.extend(vec)
    return cpp_vec


def convert_to_cpp_mat(mat):
    cpp_mat = hm.double_mat()
    for vec in mat:
        cpp_mat.append(convert_to_cpp_vec(vec))
    return cpp_mat


def to_cpp_array(array):
    dim = 0
    if type(array) in (list, np.ndarray):
        dim = 1
        if type(array[0]) in (list, np.ndarray):
            dim = 2
    if dim is 1:
        array = convert_to_cpp_vec(array)
    elif dim is 2:
        array = convert_to_cpp_mat(array)
    return array


def to_numeric(var):
    if isinstance(var, list):
        return [to_numeric(v) for v in var]
    if var == 'x' or var == sym.Symbol('x'):
        return 0
    elif var == 'y' or var == sym.Symbol('y'):
        return 1
    elif var == 'z' or var == sym.Symbol('z'):
        return 2
    else:
        return var


@convert(to_cpp_array, 'nodes', 'translation', 'dilation')
def discretize(function, nodes, translation, dilation):
    return np.array(hm.discretize(function, nodes, translation, dilation))


@convert(to_cpp_array, 'fgrid', 'nodes', 'weights')
def integrate(fgrid, nodes, weights):
    return hm.integrate(fgrid, nodes, weights)


@convert(to_cpp_array, 'fgrid', 'nodes', 'weights')
def transform(degree, fgrid, nodes, weights, forward):
    return np.array(hm.transform(degree, fgrid, nodes, weights, forward))


def triple_products(degree):
    return np.array(hm.triple_products(degree))


@convert(to_cpp_array, 'fgrid', 'nodes', 'weights')
def varf(degree, fgrid, nodes, weights):
    return np.array(hm.varf(degree, fgrid, nodes, weights))


@convert(to_cpp_array, 'var')
def varfd(dim, degree, direction, var):
    return np.array(hm.varfd(dim, degree, direction, var))


@convert(to_cpp_array, 'inp')
def tensorize(inp, dim, direction):
    return np.array(hm.tensorize(inp, dim, direction))


@convert(to_numeric, 'direction')
def _project(inp, dim, direction):
    return np.array(hm.project(inp, dim, direction))


@convert(to_cpp_array, 'inp')
def project(inp, dim, direction):
    return _project(inp, dim, direction)


def stringify(function):
    if isinstance(function, int) or isinstance(function, float):
        return str(function)
    if isinstance(function, sym.Expr):
        function = sym.ccode(function)
    if isinstance(function, str):
        function = re.sub(r'\bx\b', 'v[0]', function)
        function = re.sub(r'\by\b', 'v[1]', function)
        function = re.sub(r'\bz\b', 'v[2]', function)
        function = re.sub(r'(?<=[v])([0-9]+)', r'[\1]', function)
    return function


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


#  TODO: Ensure consistency (urbain, Thu 12 Apr 2018 09:16:39 AM BST)
def multi_indices(dim, deg_max, deg_min=0):
    return [m for m in itertools.product(range(deg_max+1), repeat=dim) if
            sum(m) <= deg_max and sum(m) >= deg_min]


def split_operator(op, func, order):
    variables = func.args
    result, rem, order = [], op.expand(), 2
    for m in multi_indices(len(variables), order):
        if rem == 0:
            result.append(0)
            continue
        test, der = 1, func
        for i, v in zip(m, variables):
            test *= v**i/math.factorial(i)
            der = sym.diff(der, v, i)
        remargs = rem.args if isinstance(rem, sym.add.Add) else [rem]
        term, rem = 0, 0
        for arg in remargs:  # Convoluted to avoid rounding errors
            termarg = arg.subs(func, test).doit()
            if termarg == 0:
                rem += arg
            else:
                term += termarg
        if isinstance(term, tuple(sym.core.all_classes)):
            term = sym.simplify(term)
        result.append(term)
    assert rem == 0
    return result


class Series:

    @staticmethod
    def natural_bissect(func, x1=0, x2=1000):
        f1, f2 = func(x1), func(x2)
        if f1 is 0:
            return x1
        elif f2 is 0:
            return x2
        assert f1*f2 < 0
        x3 = (x1+x2)//2
        f3 = func(x3)
        replace_arg = 'x2' if f1*f3 <= 0 else 'x1'
        new_args = {'x1': x1, 'x2': x2}
        new_args[replace_arg] = x3
        return Series.natural_bissect(func, **new_args)

    def __init__(self, coeffs, dim=1, mean=None, cov=None,
                 degree=None, norm=False):
        self.coeffs = coeffs/la.norm(coeffs, 2) if norm else coeffs

        self.dim = dim
        self.mean = np.zeros(self.dim) if mean is None \
            else np.asarray(mean, float)
        self.cov = np.eye(self.dim) if cov is None \
            else np.asarray(cov, float)

        eigval, eigvec = la.eig(self.cov)
        self.factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

        if degree is None:
            def obj(x):
                return int(binom(x + self.dim, x)) - len(self.coeffs)
            self.degree = Series.natural_bissect(obj)
        else:
            self.degree = degree

    @convert(to_numeric, 'direction')
    def project(self, direction):
        p_coeffs = project(self.coeffs, self.dim, direction)
        return Series(p_coeffs,
                      mean=[self.mean[direction]],
                      cov=[[self.cov[direction][direction]]])


class Quad:
    def __init__(self, nodes, weights, mean=None, cov=None):
        self.nodes = nodes
        self.weights = weights

        self.dim = len(self.nodes)
        self.mean = np.zeros(self.dim) if mean is None \
            else np.asarray(mean, float)
        self.cov = np.eye(self.dim) if cov is None \
            else np.asarray(cov, float)

        eigval, eigvec = la.eig(self.cov)
        self.factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

    @classmethod
    def gauss_hermite(cls, n_points, dim=None, mean=None, cov=None):
        if dim is not None:
            n_points = np.full(dim, n_points)
        elif isinstance(n_points, int):
            n_points = [n_points]
        nodes, weights = hermegauss_nd(n_points)
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
        function = stringify(f)
        if isinstance(function, str):
            function = discretize(function, self.nodes,
                                  self.mean, self.factor)
        return function

    def integrate(self, function):
        f_grid = self.discretize(function)
        return integrate(f_grid, self.nodes, self.weights)

    def transform(self, function, degree, norm=False):
        f_grid = self.discretize(function)
        coeffs = transform(degree, f_grid, self.nodes,
                           self.weights, forward=True)
        return Series(coeffs, self.dim, self.mean, self.cov,
                      norm=norm, degree=degree)

    def eval(self, series):
        if type(series) is np.ndarray:
            series = Series(series, self.dim, self.mean, self.cov)
        degree, coeffs = series.degree, series.coeffs
        inv = la.inv(series.factor)
        translation = inv.dot(self.mean - series.mean)
        factor = inv * self.factor
        if la.norm(factor - np.diag(np.diag(factor)), 2) > 1e-8:
            raise ValueError("Incompatible covariance matrices")
        mapped_nodes = self.nodes.copy()
        for i in range(len(self.nodes)):
            mapped_nodes[i] = self.nodes[i] * factor[i][i] + translation[i]
        return transform(degree, coeffs, mapped_nodes,
                         self.weights, forward=False)

    def varf(self, function, degree):
        f_grid = self.discretize(function)
        return varf(degree, f_grid, self.nodes, self.weights)

    @convert(to_numeric, 'directions')
    def varfd(self, function, degree, directions):
        var = self.varf(function, degree)
        eigval, _ = la.eig(self.cov)
        for dir in directions:
            var = varfd(self.dim, degree, dir, var)
            var = var/np.sqrt(eigval[dir])
        return var

    #  TODO: Improvement: tensorize when possible
    def discretize_op(self, op, func, degree, order):
        npolys = int(binom(degree + self.dim, degree))
        mat_operator = np.zeros((npolys, npolys))
        mult = list(multi_indices(self.dim, order))
        splitop = split_operator(op, func, order)
        for m, coeff in zip(mult, splitop):
            mat_operator += self.varfd(coeff, degree, ['x']*m[0] + ['y']*m[1])
        return mat_operator

    def plot(self, series, degree, factor, ax=None):
        factor = self.discretize(factor)
        if la.norm(self.cov - np.diag(np.diag(self.cov)), 2) > 1e-10:
            raise ValueError("Covariance matrix must be diagonal!")
        n_nodes = []
        r_nodes = []
        for i in range(self.dim):
            coord = 'v[{}]'.format(i)
            n_nodes.append(len(self.nodes[i]))
            r_nodes.append(self.project(i).discretize('x'))
        solution = self.eval(series)*factor
        solution = solution.reshape(*n_nodes).T
        if self.dim == 1:
            return ax.plot(*r_nodes, solution)
        elif self.dim == 2:
            return ax.contourf(*r_nodes, solution, 100)

    def factor_mapping(self):
        var = [sym.symbols('v'+ str(i)) for i in range(self.dim)]
        inv_cov = la.inv(self.cov)
        potential = 0
        for i in range(self.dim):
            for j in range(self.dim):
                potential += 0.5 * inv_cov[i][j] \
                              * (var[i] - self.mean[i]) \
                              * (var[j] - self.mean[j])
        return sym.exp(-potential/2)

    @convert(to_numeric, 'direction')
    def project(self, direction):
        return Quad([self.nodes[direction]],
                    [self.weights[direction]],
                    mean=[self.mean[direction]],
                    cov=[[self.cov[direction][direction]]])

    def series(self, coeffs, norm=False):
        return Series(coeffs,
                      dim=self.dim,
                      mean=self.mean,
                      cov=self.cov,
                      norm=norm)


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
