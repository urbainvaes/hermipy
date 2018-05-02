# TODO: Add support for sparse matrices
# TODO: Implement composite quadrature
# TODO: Ensure directions match
# TODO: Ensure variables v[0] and x can be used interchangeably
# TODO: Implement project to two dimensions
# TODO: Add function class?
# TODO: Can varfd be tensorized?

# TODO: Improve separability: for example f(x,y) * g(z)

from .cpp import hermite_cpp as hm
from scipy.special import binom
from functools import wraps
import hashlib
import math
import numpy as np
import numpy.linalg as la
import numpy.polynomial.hermite_e as herm
import re
import sympy as sym
import time

settings = {
        'cache': False,
        'cachedir': 'cache',
        'tensorize': True
        }

stats = {}


def log_stats(function):
    @wraps(function)
    def wrapper(*args, **kwargs):
        time_start = time.time()
        result = function(*args, **kwargs)
        time_end = time.time()

        key = function.__name__
        if key not in stats:
            stats[key] = {'Calls': 0, 'Time': 0}
        stats[key]['Calls'] += 1
        stats[key]['Time'] += time_end - time_start
        return result
    return wrapper


def cache(function):

    def my_hash(argument):
        if isinstance(argument, str):
            encoded_str = argument.encode('utf-8')
            return hashlib.md5(encoded_str).hexdigest()
        elif isinstance(argument, np.ndarray):
            return hashlib.md5(argument).hexdigest()
        if isinstance(argument, tuple(sym.core.all_classes)):
            return my_hash(str(argument))
        elif isinstance(argument, list):
            hashes = [my_hash(e) for e in argument]
            return my_hash(hash(frozenset(hashes)))
        elif isinstance(argument, dict):
            hashes = {kw: my_hash(argument[kw]) for kw in argument}
            return my_hash(hash(frozenset(argument)))
        elif isinstance(argument, (int, float)):
            return my_hash(str(hash(argument)))
        elif isinstance(argument, Quad):
            return my_hash(str(hash(argument)))
        elif argument is None:
            return my_hash("")
        else:
            raise ValueError("Argument type not supported")

    def error(u, v):
        if isinstance(u, (float, int)):
            return abs(u - v)
        elif isinstance(u, np.ndarray):
            return la.norm(u - v, 2)
        elif isinstance(u, tuple):
            return sum([error(ui, vi) for ui, vi in zip(u, v)])
        else:
            raise ValueError("Invalid types")

    @wraps(function)
    def wrapper(*args, **kwargs):
        use_cache = settings['cache']
        if 'cache' in kwargs:
            use_cache = kwargs['cache']
            del kwargs['cache']
        hashes, prefix = [], function.__name__
        for arg in args:
            hashes.append(my_hash(arg))
        for kw in kwargs:
            hashes.append(my_hash(kw))
            hashes.append(my_hash(kwargs[kw]))
        hash_args = my_hash(('-'.join(hashes)))

        cachedir = settings['cachedir']
        savefile = cachedir + '/' + prefix + '-' + str(hash_args) + '.npy'
        try:
            result_cache = np.load(savefile)
        except IOError:
            result = function(*args, **kwargs)
            np.save(savefile, result)
            return result

        if use_cache:
            return result_cache
        else:
            result = function(*args, **kwargs)
            assert error(result, result_cache) < 1e-10
            return result

    return wrapper


def convert_to_cpp_vec(vec):
    cpp_vec = hm.double_vec()
    cpp_vec.extend(vec)
    return cpp_vec


def convert_to_cpp_mat(mat):
    cpp_mat = hm.double_mat()
    for vec in mat:
        cpp_mat.append(convert_to_cpp_vec(vec))
    return cpp_mat


def convert_to_cpp_cube(cube):
    cpp_cube = hm.double_cube()
    for mat in cube:
        cpp_cube.append(convert_to_cpp_mat(mat))
    return cpp_cube


def to_cpp_array(*args):
    if len(args) > 1:
        return (to_cpp_array(arg) for arg in args)
    array, dim = args[0], 0
    if type(array) in (list, np.ndarray):
        dim = 1
        if type(array[0]) in (list, np.ndarray):
            dim = 2
            if type(array[0][0]) in (list, np.ndarray):
                dim = 3
    if dim is 1:
        array = convert_to_cpp_vec(array)
    elif dim is 2:
        array = convert_to_cpp_mat(array)
    elif dim is 3:
        array = convert_to_cpp_cube(array)
    return array


def to_numeric(var):

    if isinstance(var, list):
        return [to_numeric(v) for v in var]

    var_letters = {'x': 0, 'y': 1, 'z': 2}
    var_array = {'v[' + str(i) + ']': i for i in range(10)}
    var_dict = {**var_letters, **var_array}
    if var in var_dict:
        return var_dict[var]
    else:
        return var


@cache
@log_stats
def discretize(function, nodes, translation, dilation):
    nodes, translation, dilation = to_cpp_array(nodes, translation, dilation)
    return np.array(hm.discretize(function, nodes, translation, dilation))


@cache
@log_stats
def integrate(fgrid, nodes, weights):
    fgrid, nodes, weights = to_cpp_array(fgrid, nodes, weights)
    return hm.integrate(fgrid, nodes, weights)


@cache
@log_stats
def transform(degree, fgrid, nodes, weights, forward):
    fgrid, nodes, weights = to_cpp_array(fgrid, nodes, weights)
    return np.array(hm.transform(degree, fgrid, nodes, weights, forward))


@cache
@log_stats
def triple_products(degree):
    return np.array(hm.triple_products(degree))


@cache
@log_stats
def varf(degree, fgrid, nodes, weights):
    fgrid, nodes, weights = to_cpp_array(fgrid, nodes, weights)
    return np.array(hm.varf(degree, fgrid, nodes, weights))


@cache
@log_stats
def varfd(dim, degree, direction, var):
    var = to_cpp_array(var)
    return np.array(hm.varfd(dim, degree, direction, var))


@log_stats
def tensorize(inp, dim=None, direction=None):
    is_scalar = isinstance(inp[0], (float, int))
    if is_scalar and dim is None and direction is None:
        return np.prod(inp)
    inp = to_cpp_array(inp)
    if dim is not None and direction is not None:
        args = [inp, dim, direction]
    else:
        args = [inp]
    return np.array(hm.tensorize(*args))


@log_stats
def project(inp, dim, direction):
    inp = to_cpp_array(inp)
    direction = to_numeric(direction)
    return np.array(hm.project(inp, dim, direction))


@log_stats
def multi_indices(dim, degree):
    result = hm.list_multi_indices(dim, degree)
    return np.asarray(result, dtype=int)


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

    def __add__(self, other):
        assert abs(self.dim - other.dim) < 1e-8
        assert la.norm(self.mean - other.mean, 2) < 1e-8
        assert la.norm(self.cov - other.cov, 2) < 1e-8
        new_coeffs = self.coeffs + other.coeffs
        return Series(new_coeffs, dim=self.dim, mean=self.mean, cov=self.cov)

    def project(self, direction):
        direction = to_numeric(direction)
        p_coeffs = project(self.coeffs, self.dim, direction)
        return Series(p_coeffs,
                      mean=[self.mean[direction]],
                      cov=[[self.cov[direction][direction]]])


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

        self.hash = hash(frozenset({
            self.dim,
            hash(frozenset(self.nodes.flatten())),
            hash(frozenset(self.weights.flatten())),
            hash(frozenset(self.mean.flatten())),
            hash(frozenset(self.cov.flatten()))}))

    def __hash__(self):
        return self.hash

    def tensorize_at(arg_num):
        def x_ify(function):
            if not isinstance(function, tuple(sym.core.all_classes)):
                return function
            symbols = list(function.free_symbols)
            if symbols == []:
                return function
            assert len(symbols) == 1
            return function.subs(symbols[0], sym.symbols('x'))

        def split_product(expression, symbols):
            is_mul = isinstance(expression, sym.mul.Mul)
            args = expression.args if is_mul else [expression]
            result = {str(s): sym.Rational('1') for s in symbols}
            for arg in args:
                str_symbols = [str(s) for s in arg.free_symbols]
                if len(str_symbols) == 0:
                    result['v[0]'] *= arg
                elif len(str_symbols) == 1:
                    x_ified = stringify(str_symbols[0])
                    result[x_ified] *= x_ify(arg)
                else:
                    return False
            return result

        def tensorize_arg(func):
            def wrapper(*args, **kwargs):
                do_tensorize = settings['tensorize']
                if 'tensorize' in kwargs:
                    do_tensorize = kwargs['tensorize']
                    del kwargs['tensorize']
                if not do_tensorize:
                    return func(*args, **kwargs)
                function = args[arg_num]
                is_sym = isinstance(function, tuple(sym.core.all_classes))
                if not is_sym:
                    return func(*args, **kwargs)
                quad, function = args[0], function.expand()
                if isinstance(function, sym.add.Add):
                    add_terms, results = function.args, []
                    for term in add_terms:
                        new_args = list(args).copy()
                        new_args[arg_num] = term
                        func_term = func(*new_args, **kwargs)
                        results.append(func_term)
                    return sum(results[1:], results[0])
                diag_cov = np.diag(np.diag(quad.cov))
                is_diag = la.norm(quad.cov - diag_cov, 2) < 1e-10
                if quad.dim == 1 or not is_diag:
                    return func(*args, **kwargs)
                dirs = ['v[' + str(i) + ']' for i in range(quad.dim)]
                split_term = split_product(function, dirs)
                if split_term is False:
                    return func(*args, **kwargs)
                func_dirs = []
                for d in dirs:
                    new_args = list(args).copy()
                    new_args[0] = quad.project(d)
                    new_args[arg_num] = split_term[d]
                    func_dir = func(*new_args, **kwargs)
                    func_dirs.append(func_dir)
                return tensorize(func_dirs)
            return wrapper
        return tensorize_arg

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

    @tensorize_at(1)
    def integrate(self, function, l2=False):
        f_grid = self.discretize(function)
        if l2:
            w_grid = self.discretize(self.weight())
            f_grid = f_grid / w_grid
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

    @tensorize_at(1)
    def varf(self, function, degree):
        f_grid = self.discretize(function)
        return varf(degree, f_grid, self.nodes, self.weights)

    def varfd(self, function, degree, directions):
        directions = to_numeric(directions)
        var = self.varf(function, degree)
        eigval, _ = la.eig(self.cov)
        for dir in directions:
            var = varfd(self.dim, degree, dir, var)
            var = var/np.sqrt(eigval[dir])
        return var

    @cache
    def discretize_op(self, op, func, degree, order):
        npolys = int(binom(degree + self.dim, degree))
        mat_operator = np.zeros((npolys, npolys))
        mult = list(multi_indices(self.dim, order))
        splitop = split_operator(op, func, order)
        for m, coeff in zip(mult, splitop):
            mat_operator += self.varfd(coeff, degree, ['x']*m[0] + ['y']*m[1])
        return mat_operator

    #  TODO: Ensure order is right (urbain, Tue 01 May 2018)
    def plot(self, series, degree, factor, ax=None):
        factor = self.discretize(factor)
        if la.norm(self.cov - np.diag(np.diag(self.cov)), 2) > 1e-10:
            raise ValueError("Covariance matrix must be diagonal!")
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
        var = [sym.symbols('v' + str(i)) for i in range(self.dim)]
        inv_cov = la.inv(self.cov)
        potential = 0.5 * inv_cov.dot(var - self.mean).dot(var - self.mean)
        normalization = 1/(np.sqrt((2*np.pi)**self.dim * la.det(self.cov)))
        return normalization * sym.exp(-potential)

    def project(self, direction):
        direction = to_numeric(direction)
        return Quad([self.nodes[direction]],
                    [self.weights[direction]],
                    mean=[self.mean[direction]],
                    cov=[[self.cov[direction][direction]]])

    def series(self, coeffs, degree=None, norm=False):
        return Series(coeffs,
                      dim=self.dim,
                      mean=self.mean,
                      cov=self.cov,
                      degree=degree,
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
