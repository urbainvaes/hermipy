## @Package spectral
# Hermite spectral method
#
# This module provides useful functions for the Hermite spectral method.
# Here we deal exclusively with normalized Hermite polynomials, in the
# sense that:
# \f[
#     \int_{\mathbb R^n} H^\alpha(x) \, H^\beta(x) \, \rho(x) \, \mathrm d x = 1
# \f]
import collections
import itertools
# import numba as nb
import numpy as np
import numpy.linalg as la
import numpy.polynomial.hermite_e as herm
import scipy.special

from libhermite import hermite as hm


def convert_to_cpp_vec(vec):
    cpp_vec = hm.double_vec()
    cpp_vec.extend(vec)
    return cpp_vec


def convert_to_cpp_mat(mat):
    cpp_mat = hm.double_mat()
    for vec in mat:
        cpp_mat.append(convert_to_cpp_vec(vec))
    return cpp_mat


def convert_to_numpy_vec(vec):
    numpy_vec = np.zeros(len(vec))
    for i in range(len(vec)):
        numpy_vec[i] = vec[i]
    return numpy_vec


def convert_to_numpy_mat(mat):
    numpy_mat = np.zeros((len(mat), len(mat[0])))
    for i in range(len(mat)):
        numpy_mat[i] = convert_to_numpy_vec(mat[i])
    return numpy_mat


def convert_to_list(array):
    if not isinstance(array, collections.Iterable):
        return array
    else:
        result = []
        for elem in result:
            result.append(convert_to_list(elem))


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
def hermegauss_nd(deg, dim=None):
    if dim is not None:
        n_nodes = np.full(dim, deg)
    elif isinstance(deg, int):
        n_nodes = [deg]
    else:
        n_nodes = deg
    dim = len(n_nodes)
    nodes_multidim = []
    weights_multidim = []
    for i in range(dim):
        nodes_1d, weights_1d = herm.hermegauss(n_nodes[i])
        weights_1d = weights_1d/np.sqrt(2*np.pi)  # Normalize
        nodes_multidim.append(nodes_1d)
        weights_multidim.append(weights_1d)
    return nodes_multidim, weights_multidim

## Compute the weighted multi-dimensional integral of a function using a
# Gauss-Hermite quadrature
#
# See spectral::hermegauss_nd for a full documentation of the parameters
#
# ## Note
# It is important that this function be able to take a list of
def integrate_simple_quad(f, nodes, weights, mean=None, cov=None):

    dim = len(nodes)
    if mean is None:
        mean = np.zeros(dim)
    if cov is None:
        cov = np.eye(dim)

    eigval, eigvec = la.eig(cov)
    mat_factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

    cpp_translation = convert_to_cpp_vec(mean)
    cpp_dilation = convert_to_cpp_mat(mat_factor)

    cpp_nodes = convert_to_cpp_mat(nodes)
    cpp_weights = convert_to_cpp_mat(weights)

    if isinstance(f, str):
        f_nodes = hm.discretize(f, cpp_nodes, cpp_translation, cpp_dilation)
    else:
        f_nodes = f

    result = hm.integrate(f_nodes, cpp_nodes, cpp_weights)
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
def integrate(f, deg, dim=None, mean=None, cov=None):
    nodes, weights = hermegauss_nd(deg, dim)
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
def hermval(x, c, mean=None, cov=None):
    dim = len(x[0])
    if mean is None:
        mean = np.zeros(dim)
    if cov is None:
        cov = np.eye(dim)
    eigval, eigvec = la.eig(cov)
    mat_factor = np.matmul(eigvec,np.sqrt(np.diag(eigval)))
    x_standard = la.inv(mat_factor)(x - mean).T


def monval(x, c, mean=None, cov=None):
    dim = len(x[0])
    if mean is None:
        mean = np.zeros(dim)
    if cov is None:
        cov = np.eye(dim)
    eigval, eigvec = la.eig(cov)
    mat_factor = np.matmul(eigvec,np.sqrt(np.diag(eigval)))
    x_standard = la.inv(mat_factor)(x - mean).T


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

def transform_simple_quad(f, degree, nodes, weights, mean=None, cov=None):

    dim = len(nodes)
    if mean is None:
        mean = np.zeros(dim)
    if cov is None:
        cov = np.eye(dim)

    eigval, eigvec = la.eig(cov)
    mat_factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

    cpp_translation = convert_to_cpp_vec(mean)
    cpp_dilation = convert_to_cpp_mat(mat_factor)

    n_poly = int(scipy.special.binom(degree + dim, dim))
    result = np.zeros(n_poly)

    cpp_nodes = convert_to_cpp_mat(nodes)
    cpp_weights = convert_to_cpp_mat(weights)

    if isinstance(f, str):
        f_nodes = hm.discretize(f, cpp_nodes, cpp_translation, cpp_dilation)
    else:
        f_nodes = f

    result = hm.transform(degree, f_nodes, cpp_nodes, cpp_weights, True)
    return convert_to_numpy_vec(result)


def transform_composite_quad(f, degree, nodes, weights, mean=None, cov=None):

    dim = len(nodes[0])
    n_poly = int(scipy.special.binom(degree + dim, dim))
    result = np.zeros(n_poly)

    for i in range(len(nodes)):
        f_arg = f if isinstance(f, str) else f[i]
        result += transform_simple_quad(f_arg, degree, nodes[i],
                                        weights[i], mean, cov)
    return result


## Enumerate the N-dimensional multi-indices \f$ \alpha \f$ such that \f$ |\alpha| \leq b \f$.
#
# \param dim dimension of the multi-indices
# \param upper upper bound of multi-indices
def multi_indices(dim, deg_max, deg_min=0):
    return [m for m in itertools.product(range(deg_max+1), repeat=dim) if \
            sum(m) <= deg_max and sum(m) >= deg_min]

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
    def __init__(self, deg, dim=None, mean=None, cov=None):

        # Nodes and weights of the quadrature
        nodes, weights = hermegauss_nd(deg, dim)

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
    def integrate(self, f, mean=None, cov=None):
        return integrate_composite_quad(f, self.nodes, self.weights, self.mean, self.cov)
