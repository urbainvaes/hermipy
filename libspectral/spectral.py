## @Package spectral
# Hermite spectral method
#
# Provide useful functions for the finite element method
import numpy as np
import numpy.linalg as la
import numpy.polynomial.hermite_e as herm


## Compute an N-dimensional Gauss-hermite quadrature
#
# This function returns the sample points and weights of a Gauss-Hermite
# quadrature with the following general weight
# \f[
# g(\mu;\Sigma, x) := \frac{1}{{\sqrt{(2\pi)^k|\boldsymbol\Sigma|}}} \exp\left(-\frac 1 2 ({\mathbf x}-{\boldsymbol\mu})^\mathrm{T}{\boldsymbol\Sigma}^{-1}({\mathbf x}-{\boldsymbol\mu})\right),
# \f]
# where \f$ \mu \f$ = \p mean and \f$ \Sigma \f$ = cov.
# The number of sample points to use along each eigen-direction of covariance
# matrix is calculated from \p deg and \p dim as follows:
# \code{.py}
# if dim is not None:
#     n_nodes = np.full(dim, deg)
# elif isinstance(deg, int):
#     n_nodes = [deg]
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

def hermegauss_nd(deg, dim=None, mean=None, cov=None):
    if dim is not None:
        n_nodes = np.full(dim, deg)
    elif isinstance(deg, int):
        n_nodes = [deg]
    dim = len(n_nodes)
    if mean is None:
        mean = np.zeros(dim)
    if cov is None:
        cov = np.eye(dim)
    nodes_multidim = []
    weights_multidim = []
    for i in range(dim):
        nodes_1d, weights_1d = herm.hermegauss(n_nodes[i])
        weights_1d = weights_1d/np.sqrt(2*np.pi)  # Normalize
        nodes_multidim.append(nodes_1d)
        weights_multidim.append(weights_1d)
    grid_nodes = np.meshgrid(*nodes_multidim)
    grid_weights = np.meshgrid(*weights_multidim)
    nodes_flattened = []
    weights_flattened = []
    for i in range(dim):
        nodes_flattened.append(grid_nodes[i].flatten())
        weights_flattened.append(grid_weights[i].flatten())
    nodes_standard = np.vstack(nodes_flattened).T
    weights = np.prod(np.vstack(weights_flattened), axis=0)
    eigval, eigvec = la.eig(cov)
    mat_factor = np.matmul(eigvec,np.sqrt(np.diag(eigval)))
    nodes = np.matmul(mat_factor, nodes_standard.T).T + mean
    return nodes, weights


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

# \param f The function to integrate. # \returns the value of the multi-dimensional integral
def integrate(f, deg, dim=None, mean=None, cov=None, normalize=True):
    nodes, weights = hermegauss_nd(deg, dim, mean, cov, normalize)
    return np.dot(f(nodes.T), weights)


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
    def __init__(self, deg, dim=None, mean=None, cov=None, normalize=True):

        # Nodes and weights of the quadrature
        nodes, weights = hermegauss_nd(deg, dim, mean, cov, normalize)

        ## The sample points of the quadrature.
        self.nodes = nodes

        ## The weights of the quadrature.
        self.weights = weights

    ## Compute a multi-dimensional integral
    #
    # See spectral::integrate for a description of the parameters.
    def integrate(self, f):
        return np.dot(f(self.nodes.T), self.weights)
