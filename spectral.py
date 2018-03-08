## @Package spectral
# Hermite spectral method
#
# This module provides useful functions for the Hermite spectral method.
# Here we deal exclusively with normalized Hermite polynomials, in the
# sense that:
# \f[
#     \int_{\mathbb R^n} H^\alpha(x) \, H^\beta(x) \, \rho(x) \, \mathrm d x = 1
# \f]
# Import modules {{{
import collections
import itertools
import numpy as np
import numpy.linalg as la
import numpy.polynomial.hermite_e as herm
import scipy.special
import sympy as sy
import re

from libhermite import hermite as hm
from inspect import signature
# }}}
# Functions for compatibility with C++ {{{
def convert_to_cpp_vec(vec):
    cpp_vec = hm.double_vec()
    cpp_vec.extend(vec)
    return cpp_vec


def convert_to_cpp_mat(mat):
    cpp_mat = hm.double_mat()
    for vec in mat:
        cpp_mat.append(convert_to_cpp_vec(vec))
    return cpp_mat


def convert_to_cpp_array(array):
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

# }}}
# Define decorators to ease information exchange with C++ {{{
def args_to_cpp(*names):
    def convert(function):
        sig = signature(function)

        def wrapper(*args, **kwargs):
            ba = sig.bind(*args, **kwargs)
            args_od = ba.arguments
            for key in args_od:
                if key in names:
                    arg = args_od[key]
                    args_od[key] = convert_to_cpp_array(arg)
            return function(**args_od)
        return wrapper
    return convert
# }}}
## Compute an N-dimensional Gauss-hermite quadrature
#
# This function returns the sample points and weights of a Gauss-Hermite
# quadrature with the following general weight
# \f[
#     g(\mu;\Sigma, x) := \frac{1}{{\sqrt{(2\pi)^k|\boldsymbol\Sigma|}}} \exp\left(-\frac 1 2 ({\mathbf x}-{\boldsymbol\mu})^\mathrm{T}{\boldsymbol\Sigma}^{-1}({\mathbf x}-{\boldsymbol\mu})\right),
# \f]
# where \f$ \mu \f$ = \p mean and \f$ \Sigma \f$ = cov.
# The number of sample points to use along each eigen-direction of covariance
# matrix is calculated from \p deg and \p dim as follows:
# \code{.py}
# if dim is not None:
#     n_nodes = np.full(dim, deg)
# elif isinstance(deg, int):
#     n_nodes = [deg]
# else:
#     n_nodes = deg
# \endcode
#
# ## Usage
# In most casse, one will want quadratures that have the same number of sample
# points along each eigen-direction. For sparse grid methods, however, it is
# necessary to construct quadratures that do not.
#
# ## Precision
# The weights and sample points returned can be used to calculate exacly the
# weighted integral of polynomials of degree up to 2*\p n_nodes - 1.
# \param deg see above.
# \param dim see above.
# \param mean the mean of the Gaussian weight.
# \param cov the covariance matrix of the Gaussian weight.
# \retval nodes A 2-D array containing the sample points of the quadrature,
# with `nodes[0]` containing the coordinates of the first node, etc.
# \retval weights  A 1-D array containing the weights of the quadrature.
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


def stringify(function): # {{{
    if not isinstance(function, str):
        function = sy.ccode(function)

    # If x, y, z notation is used
    subst_vars = function
    subst_vars = re.sub(r'\bx\b', 'v[0]', subst_vars)
    subst_vars = re.sub(r'\by\b', 'v[1]', subst_vars)
    subst_vars = re.sub(r'\bz\b', 'v[2]', subst_vars)

    # If v0, v1, ... notation is used
    subst_vars = re.sub(r'(?<=[v])([0-9]+)', r'[\1]', subst_vars)

    return subst_vars # }}}


## Discretize a function on a simple quadrature
@args_to_cpp('nodes')
def discretize(function, nodes, mean=None, cov=None):

    is_iterable = isinstance(function, collections.Iterable)
    is_string = isinstance(function, str)

    if is_iterable and not is_string:
        return function

    function = stringify(function)

    dim = len(nodes)
    if mean is None:
        mean = np.zeros(dim)
    if cov is None:
        cov = np.eye(dim)

    eigval, eigvec = la.eig(cov)
    mat_factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

    cpp_translation = convert_to_cpp_vec(mean)
    cpp_dilation = convert_to_cpp_mat(mat_factor)
    cpp_discretization = hm.discretize(function, nodes,
                                       cpp_translation, cpp_dilation)
    return np.array(cpp_discretization)


## Compute the weighted multi-dimensional integral of a function using a
# Gauss-Hermite quadrature
#
# See spectral::hermegauss_nd for a full documentation of the parameters
#
# ## Note
# It is important that this function be able to take a list of
@args_to_cpp('nodes', 'weights')
def integrate_simple_quad(f, nodes, weights, mean=None, cov=None):
    f_grid = convert_to_cpp_vec(discretize(f, nodes, mean, cov))
    result = hm.integrate(f_grid, nodes, weights)
    return result


def integrate_composite_quad(f, nodes, weights, mean=None, cov=None):
    result = 0
    for i in range(len(nodes)):
        f_arg = f if isinstance(f, str) else f[i]
        result += integrate_simple_quad(f_arg, nodes[i],
                                        weights[i], mean, cov)
    return result


## Compute the weighted multi-dimensional integral of a function using a
# Gauss-Hermite quadrature
#
# See spectral::hermegauss_nd for a full documentation of the parameters
#
# \param f The function to integrate. # \returns the value of the multi-dimensional integral
def integrate(f, n_points, dim=None, mean=None, cov=None):
    nodes, weights = hermegauss_nd(n_points, dim)
    return integrate_simple_quad(f, nodes, weights, mean=mean, cov=cov)


## Convert a hermite series to a polynomial
#
# This function is the same as the one provided in the module
# `numpy.polynomial.hermite_e`, except that the input coefficients correspond
# to hermite polynomials normalized to have norm one with respect to the
# standard normal weight.
def herm_to_poly(c):
    herme_coeffs = c/np.sqrt(np.sqrt(2*np.pi)*np.arange(len(c)))
    return herm.herme2poly(herme_coeffs)


## Evaluate a normalized Hermite series at points x.
# \param c an array of coefficients ordered so that c[n]
# contains the coefficient of degree n.
# If c is multidimensional the remaining indices enumerate
# multiple polynomials. In the two dimensional case the
# coefficients may be thought of as stored in the columns of
# c.

## Expand Hermite function with complexity O(n*n)
#
# This algorithm is based on the following recursion formula for
# probabilist Hermite polynomial:
# \f[
#     He_{n+1}(x) = x \, He_{n} (x) - n \, He_{n-1}.
# \f]
# Since in one dimension the normalized Hermite polynomials are related to the
# probabilist Hermite polynomials by \f$ H_n = He_n / \sqrt{n!} \f$, the
# recursion formula becomes
# \f[
#     \sqrt{n+1} \, H_{n+1}(x) = x \, H_n(x) - \sqrt{n} H_{n-1}(x).
# \f]
# This allows the evaluation of the \f$ n \f$ first Hermite polynomials
# in one point in \f$ \mathcal O(n) \f$ operations, compared to \f$
# \mathcal O(n^2) \f$ for a naive implementation.

@args_to_cpp('nodes', 'weights')
def transform_simple_quad(f, degree, nodes, weights, mean=None,
                          cov=None, forward=True):

    f_grid = convert_to_cpp_vec(discretize(f, nodes, mean, cov))
    dim = len(nodes)
    n_poly = int(scipy.special.binom(degree + dim, dim))
    result = np.zeros(n_poly)

    result = hm.transform(degree, f_grid, nodes, weights, forward)
    return np.array(result)


def transform_composite_quad(f, degree, nodes, weights,
                             mean=None, cov=None, forward=True):

    dim = len(nodes[0])
    n_poly = int(scipy.special.binom(degree + dim, dim))
    result = np.zeros(n_poly)

    for i in range(len(nodes)):
        #  FIXME: desc (urbain, Wed 28 Feb 2018 02:47:04 PM GMT)
        # f_arg = f if isinstance(f, str) else f[i]
        result += transform_simple_quad(f, degree, nodes[i],
                                        weights[i], mean, cov, forward)
    return result


def eval_simple_quad(coeffs, degree, nodes):
    cpp_coeffs = convert_to_cpp_vec(coeffs)
    return transform_simple_quad(cpp_coeffs, degree, nodes, nodes, forward=False)


def eval_composite_quad(coeffs, degree, nodes):
    result = []

    for i in range(len(nodes)):
        result_quad = eval_simple_quad(coeffs, degree, nodes[i])
        result.append(result_quad)
    return result

## Enumerate the N-dimensional multi-indices \f$ \alpha \f$ such that \f$ |\alpha| \leq b \f$.
#
# \param dim dimension of the multi-indices
# \param upper upper bound of multi-indices
def multi_indices(dim, deg_max, deg_min=0):
    return [m for m in itertools.product(range(deg_max+1), repeat=dim) if \
            sum(m) <= deg_max and sum(m) >= deg_min]


@args_to_cpp('nodes', 'weights')
def varf_simple_quad(f, degree, nodes, weights, mean, cov):
    f_grid = convert_to_cpp_vec(discretize(f, nodes, mean, cov))
    result = hm.varf(degree, f_grid, nodes, weights)
    return np.array(result)

## Gauss-Hermite quadrature.
#
# To compute several integrals using the same quadrature, it
# is advantageous to create a `Quad` object and use the class
# function `integrate`. This way, the quadrature needs not be
# recalculated at every integration.
class Quad:
    ## A Gauss-Hermite quadrature
    #
    # See spectral::hermegauss_nd for a description of the parameters.
    def __init__(self, n_points, dim=None, mean=None, cov=None):

        # Nodes and weights of the quadrature
        nodes, weights = hermegauss_nd(n_points, dim)

        ## The sample points of the quadrature.
        self.nodes = [nodes]

        ## The weights of the quadrature.
        self.weights = [weights]

        ## Mean of the Gaussian weight
        self.mean = mean

        ## Covariance of the Gaussian weight
        self.cov = cov


    ## Compute a multi-dimensional integral
    #
    # See spectral::integrate for a description of the parameters.
    def integrate(self, f):
        return integrate_composite_quad(f, self.nodes, self.weights, self.mean, self.cov)

    def discretize(self, f):
        return discretize(f, self.nodes[0], self.mean, self.cov)

    def transform(self, function, degree):
        return transform_composite_quad(function, degree, self.nodes,
                                        self.weights, self.mean, self.cov)

    def eval(self, coeffs, degree):
        return eval_composite_quad(coeffs, degree, self.nodes)

    def varf(self, function, degree):
        return varf_simple_quad(function, degree, self.nodes[0],
                                self.weights[0], self.mean, self.cov)


class hermite_series:

    def __init__(self, coeffs, mean=None, cov=None):
        self.coeffs = coeffs
        self.mean = mean
        self.cov = cov

    def eval(self, degree, nodes):
        return eval_simple_quad(self.coeffs, degree, nodes)
