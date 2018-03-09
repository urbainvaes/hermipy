from libhermite import hermite_python as hm

import numpy as np
import numpy.linalg as la
import numpy.polynomial.hermite_e as herm
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


class Series:

    def __init__(self, coeffs, dim=1, mean=None, cov=None):
        self.coeffs = coeffs
        self.mean = np.zeros(dim) if mean is None else mean
        self.cov = np.eye(dim) if cov is None else cov

        eigval, eigvec = la.eig(self.cov)
        self.factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))


class Quad:
    def __init__(self, nodes, weights, mean=None, cov=None):
        self.nodes = nodes
        self.weights = weights

        self.dim = len(self.nodes)
        self.mean = np.zeros(self.dim) if mean is None else np.array(mean)
        self.cov = np.eye(self.dim) if cov is None else np.array(cov)

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

    #  TODO: Add mapping (urbain, Thu 08 Mar 2018 11:07:49 PM GMT)
    def mapped_nodes(self):
        grid_nodes = np.meshgrid(self.nodes)
        nodes_flattened = []
        for i in range(len(self.nodes)):
            nodes_flattened.append(grid_nodes[i].flatten())
        nodes = np.vstack(nodes_flattened)
        return nodes

    def discretize(self, f):
        function = stringify(f)
        if isinstance(function, str):
            function = hm.discretize(function, self.nodes,
                                     self.mean, self.factor)
        return function

    def integrate(self, function):
        f_grid = self.discretize(function)
        return hm.integrate(f_grid, self.nodes, self.weights)

    def transform(self, function, degree):
        f_grid = self.discretize(function)
        coeffs = hm.transform(degree, f_grid, self.nodes,
                              self.weights, forward=True)
        return Series(coeffs, self.dim, self.mean, self.cov)

    def eval(self, series, degree):
        if type(series) is np.ndarray:
            series = Series(series, self.dim, self.mean, self.cov)
        coeffs = series.coeffs
        translation = self.mean - series.mean
        factor = la.inv(series.factor) * self.factor
        if la.norm(factor - np.diag(np.diag(factor)), 2) > 1e-8:
            raise ValueError("Incompatible covariance matrices")
        mapped_nodes = self.nodes.copy()
        for i in range(len(self.nodes)):
            mapped_nodes[i] = self.nodes[i]*factor[i] + translation[i]
        return hm.transform(degree, coeffs, mapped_nodes,
                            self.weights, forward=False)

    def varf(self, function, degree):
        f_grid = self.discretize(function)
        return hm.varf(degree, f_grid, self.nodes, self.weights)


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
