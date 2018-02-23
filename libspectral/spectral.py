## @Package spectral
# Hermite spectral method
#
# This module provides useful functions for the Hermite spectral method.
# Here we deal exclusively with normalized Hermite polynomials, in the
# sense that:
# \f[
#     \int_{\mathbb R^n} H^\alpha(x) \, H^\beta(x) \, \rho(x) \, \mathrm d x = 1
# \f]
import itertools
import numba as nb
import numpy as np
import numpy.linalg as la
import numpy.polynomial.hermite_e as herm
import scipy.special

from libhermite import hermite as hm


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
    return [nodes_multidim], [weights_multidim]


## Compute the weighted multi-dimensional integral of a function using a
# Gauss-Hermite quadrature
#
# See spectral::hermegauss_nd for a full documentation of the parameters
#
# ## Note
# It is important that this function be able to take a list of
def integrate_with_quad(f, nodes, weights, mean=None, cov=None):

    dim = len(nodes[0])
    if mean is None:
        mean = np.zeros(dim)
    if cov is None:
        cov = np.eye(dim)

    eigval, eigvec = la.eig(cov)
    mat_factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

    cpp_nodes = hm.double_cube()
    for mat in nodes:
        m = hm.double_mat()
        for vec in mat:
            v = hm.double_vec()
            v.extend(vec)
            m.append(v)
        cpp_nodes.append(m)

    cpp_weights = hm.double_cube()
    for mat in weights:
        m = hm.double_mat()
        for vec in mat:
            v = hm.double_vec()
            v.extend(vec)
            m.append(v)
        cpp_weights.append(m)

    cpp_translation = hm.double_vec()
    cpp_translation.extend(mean)

    cpp_dilation = hm.double_mat()
    for vec in mat_factor:
        v = hm.double_vec()
        v.extend(vec)
        cpp_dilation.append(v)

    return hm.integrate_from_string(f, cpp_nodes, cpp_weights,
                                    cpp_translation, cpp_dilation)

## Compute the weighted multi-dimensional integral of a function using a
# Gauss-Hermite quadrature
#
# See spectral::hermegauss_nd for a full documentation of the parameters
#
# ## Note
# It is important that this function be able to take a list of
# vectors as input. For instance, the following function would
# not work
# \code{.py}
# def function(vector):
#     return 1
# \endcode
# An easy workaround in this case is to write
# \code{.py}
# def function(vector):
#     return 0*v[0] + 1
# \endcode
#
# \param f The function to integrate. # \returns the value of the multi-dimensional integral
def integrate(f, deg, dim=None, mean=None, cov=None):
    nodes, weights = hermegauss_nd(deg, dim)
    return integrate_with_quad(f, nodes, weights, mean=mean, cov=cov)


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
def herm_transform(function, dim, deg, n_points, mean=None, cov=None):

    @nb.jit(fastmath=True)
    def numba_f(v):
        return function(v)

    quad_unit = Quad(n_points, dim)
    quad_real = Quad(n_points, dim, mean, cov)
    n_poly = int(scipy.special.binom(deg + dim, dim))
    mult_inds = multi_indices(dim, deg)
    exp_coeffs = {}
    for m in mult_inds:
        exp_coeffs[m] = 0
    h = np.zeros((dim, deg+1))
    sq = np.sqrt(range(deg+1))
    for d in range(dim):
        h[d][0] = 1

    @nb.jit(nopython=True, fastmath=True)
    def herm_eval_1d(x, deg):
        result = np.zeros(deg + 1)
        result[0] = 1
        result[1] = x
        for n in range(1, deg):
            result[n+1] = (1/sq[n+1])*(x*result[n]-sq[n]*result[n-1])
            # result[n+1] = (x*result[n]-n*result[n-1])
        return result

    for i_node in range(len(quad_unit.nodes)):
        u_node = quad_unit.nodes[i_node]
        r_node = quad_real.nodes[i_node]
        for d in range(dim):
            h[d] = herm_eval_1d(u_node[d], deg)
        for m in mult_inds:
            vals = [h[d][m[d]] for d in range(dim)]
            # print(mult_inds[i_mult])
            # print(vals)
            exp_coeffs[m] += np.prod(vals)*numba_f(r_node)*quad_unit.weights[i_node]
    return exp_coeffs

def inv_hermite_transform(coeffs, dim, n_points):
    quad = Quad(n_points, dim)

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
        nodes, weights = hermegauss_nd(deg, dim, mean, cov)

        ## The sample points of the quadrature.
        self.nodes = nodes

        ## The weights of the quadrature.
        self.weights = weights

    ## Compute a multi-dimensional integral
    #
    # See spectral::integrate for a description of the parameters.
    def integrate(self, f, mean=None, cov=None):
        return integrate_with_quad(f, self.nodes, self.weights, mean=mean, cov=cov)
