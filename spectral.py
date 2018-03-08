from libhermite import hermite_python as hm
import itertools
import numpy as np
import numpy.linalg as la
import numpy.polynomial.hermite_e as herm
import scipy.special
import sympy as sy
import re


def stringify(function):
    if isinstance(function, sy.Expr):
        function = sy.ccode(function)
    if isinstance(function, str):
        function = re.sub(r'\bx\b', 'v[0]', function)
        function = re.sub(r'\by\b', 'v[1]', function)
        function = re.sub(r'\bz\b', 'v[2]', function)
        function = re.sub(r'(?<=[v])([0-9]+)', r'[\1]', function)
    return function


def hermegauss_nd(n_points, dim=None):
    if dim is not None:
        n_nodes = np.full(dim, n_points)
    elif isinstance(n_points, int):
        n_nodes = [n_points]
    else:
        n_nodes = n_points
    dim = len(n_nodes)
    nodes_multidim = []
    weights_multidim = []
    for i in range(dim):
        nodes_1d, weights_1d = herm.hermegauss(n_nodes[i])
        weights_1d = weights_1d/np.sqrt(2*np.pi)  # Normalize
        nodes_multidim.append(nodes_1d)
        weights_multidim.append(weights_1d)
    return nodes_multidim, weights_multidim


class Quad:
    def __init__(self,
                 n_points=None,
                 dim=None,
                 nodes=None,
                 weights=None,
                 mean=None,
                 cov=None):

        if n_points is not None:
            nodes, weights = hermegauss_nd(n_points, dim)

        self.nodes = nodes
        self.weights = weights

        dim = len(self.nodes)
        self.mean = np.zeros(dim) if mean is None else mean
        self.cov = np.eye(dim) if cov is None else cov

    def discretize(self, f):
        function = stringify(f)
        if isinstance(function, str):
            eigval, eigvec = la.eig(self.cov)
            factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))
            function = hm.discretize(function, self.nodes, self.mean, factor)
        return function

    def integrate(self, function):
        f_grid = self.discretize(function)
        return hm.integrate(f_grid, self.nodes, self.weights)

    def transform(self, function, degree):
        f_grid = self.discretize(function)
        return hm.transform(degree, f_grid, self.nodes, self.weights, forward=True)

    def eval(self, coeffs, degree):
        return hm.transform(degree, coeffs, self.nodes, self.weights, forward=False)

    def varf(self, function, degree):
        f_grid = self.discretize(function)
        return hm.varf(degree, f_grid, self.nodes, self.weights)


class CompQuad:
    def __init__(self, quads, weights):
        self.quads = quads
        self.weights = weights

    # def integrate(f):


class hermite_series:

    def __init__(self, coeffs, mean=None, cov=None):
        self.coeffs = coeffs
        self.mean = mean
        self.cov = cov

    def eval(self, degree, nodes):
        return eval_simple_quad(self.coeffs, degree, nodes)

# def herm_to_poly(c):
#     herme_coeffs = c/np.sqrt(np.sqrt(2*np.pi)*np.arange(len(c)))
#     return herm.herme2poly(herme_coeffs)
