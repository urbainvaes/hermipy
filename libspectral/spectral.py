## @Package spectral
# Hermite spectral method
#
# Provide useful functions for the finite element method
import numpy as np
import numpy.polynomial.hermite_e as herm

## Compute an N-dimensional Gauss-hermite quadrature
#
# \param degrees list containing the number of quadrature nodes in each dimension.
# \param normalize see `Quad::integrate`
# \returns the nodes and weights of the N-dimensional quadrature.
def hermegauss_nd(degrees, normalize=True):
    dim = len(degrees)
    nodes_multidim = []
    weights_multidim = []
    for i in range(dim):
        nodes_1d, weights_1d = herm.hermegauss(degrees[i])
        if normalize:  # Normalize
            weights_1d = weights_1d/np.sqrt(2*np.pi)
        nodes_multidim.append(nodes_1d)
        weights_multidim.append(weights_1d)
    grid_nodes = np.meshgrid(*nodes_multidim)
    grid_weights = np.meshgrid(*weights_multidim)
    nodes_flattened = []
    weights_flattened = []
    for i in range(dim):
        nodes_flattened.append(grid_nodes[i].flatten())
        weights_flattened.append(grid_weights[i].flatten())
    nodes = np.vstack(nodes_flattened)
    weights = np.prod(np.vstack(weights_flattened), axis=0)
    return nodes, weights


## Compute the multi-dimensional integral of a function
#
# \param f: the function to integrate.
# \param dimension: the dimension of the domain of `f`.
# \param degree: the number of nodes to use for the Gauss-Hermite integral.
def integrate(f, dimension, degree):
    nodes, weights = hermegauss_nd(np.full(dimension, degree))
    # f_nodes = np.apply_along_axis(f, 0, nodes)  # Sometimes, f(nodes) work
    f_nodes = f(nodes)
    return np.dot(f_nodes, weights)


## Gauss-Hermite quadrature.
#
# To compute several integrals using the same quadrature, it
# is advantageous to create a `Quad` object and use the class
# function `integrate`. This way, the quadrature needs not be
# recalculated at every integration.
class Quad:
    ## Create a `Quad` object
    #
    # \param degrees a list containing the number of quadrature nodes to use in
    # each dimension.
    # \param normalize a boolean determining whether or not to include
    # the prefactor \f$ 1/\sqrt{2\pi} \f$ in the weight function.
    # When \p normalize is True, the weight used is the normal density with
    # mean 0 and identity covariance:
    # \f[
    #      \frac{1}{(\sqrt{2\pi})^n} e^{-\frac{x^T x}{2}},
    # \f]
    # where \f$ n \f$ is the length of \p degrees.
    def __init__(self, degrees, normalize=True):
        nodes, weights = hermegauss_nd(degrees, normalize)

        ## The nodes of the quadrature.
        self.nodes = nodes

        ## The weights of the quadrature.
        self.weights = weights

    ## Compute integral
    # \param f The function to integrate. It is important that
    # this function be able to take a list of vectors as
    # input. For instance, the following function would not work
    # \code{.py}
    # def function(vector):
    #     return 1
    # \endcode
    # An easy workaround in this case is to write
    # \code{.py}
    # def function(vector):
    #     return 0*v[0] + 1
    # \endcode
    # \returns The multi-dimensional integral using the quadrature
    def integrate(self, f):
        f_nodes = f(self.nodes)
        return np.dot(f_nodes, self.weights)
