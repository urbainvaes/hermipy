## @Package spectral
# Hermite spectral method
#
# This module provides useful functions for the Hermite spectral method.
# Here we deal exclusively with normalized Hermite polynomials, in the
# sense that:
# \f[
#     \int_{\mathbb R^n} H^\alpha(x) \, H^\beta(x) \, \rho(x) \, \mathrm d x = 1
# \f]


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


## Compute the weighted multi-dimensional integral of a function using a
# Gauss-Hermite quadrature
#
# See spectral::hermegauss_nd for a full documentation of the parameters
#
# ## Note
# It is important that this function be able to take a list of

## Compute the weighted multi-dimensional integral of a function using a
# Gauss-Hermite quadrature
#
# See spectral::hermegauss_nd for a full documentation of the parameters
#
# \param f The function to integrate. # \returns the value of the multi-dimensional integral

## Convert a hermite series to a polynomial
#
# This function is the same as the one provided in the module
# `numpy.polynomial.hermite_e`, except that the input coefficients correspond
# to hermite polynomials normalized to have norm one with respect to the
# standard normal weight.

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


## Enumerate the N-dimensional multi-indices \f$ \alpha \f$ such that \f$ |\alpha| \leq b \f$.
#
# \param dim dimension of the multi-indices
# \param upper upper bound of multi-indices

## Gauss-Hermite quadrature.
#
# To compute several integrals using the same quadrature, it
# is advantageous to create a `Quad` object and use the class
# function `integrate`. This way, the quadrature needs not be
# recalculated at every integration.
