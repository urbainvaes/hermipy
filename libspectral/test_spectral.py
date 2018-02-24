import unittest
import numpy as np
import spectral as sp


class TestIntegrate(unittest.TestCase):

    ## Test that the weighted integral of the constant function
    # \f$ f := 1 \f$ is equal to 1 for different quadrature
    # sizes.
    def test_normalization_nodes(self):
        deg = [2**i for i in range(8)]
        for i in deg:
            integral = sp.integrate('1', [i])
            self.assertAlmostEqual(integral, 1)

            integral = sp.integrate('v[0]', [i])
            self.assertAlmostEqual(integral, 0)

    ## Test that the weighted integral of the constant function
    # \f$ f := 1 \f$ is equal to 1 in dimensions 1,2 and 3.
    def test_normalization_dim(self):
        for i in range(1, 4):
            integral = sp.integrate('1', [i])
            self.assertAlmostEqual(integral, 1)

            integral = sp.integrate('v[0]', 100, dim=i)
            self.assertAlmostEqual(integral, 0)

            integral = sp.integrate('v[0]*v[0]', 100, dim=i)
            self.assertAlmostEqual(integral, 1)

    ## Test that the first moment of the weight has the correct value.
    def test_mean(self):
        dim = 3
        mean = np.random.random(dim)
        for i in range(len(mean)):
            fun = 'v[{}]'.format(i)
            coord = sp.integrate(fun, 8, dim=dim, mean=mean)
            self.assertAlmostEqual(coord, mean[i])

    ## Test that the second moment of the weight has the correct value.
    def test_covariance(self):
        dim = 3
        rand_mat = np.random.random((dim, dim))
        cov = np.matmul(rand_mat.T, rand_mat)
        for i in range(len(cov)):
            for j in range(len(cov)):
                fun = 'v[{}]*v[{}]'.format(i, j)
                cov_ij = sp.integrate(fun, 8, dim=dim, cov=cov)
                self.assertAlmostEqual(cov_ij, cov[i][j])

    ## Test the quadrature object
    def test_quad(self):
        dim = 3
        rand_mat = np.random.random((dim, dim))
        mean = np.random.random(dim)
        cov = np.matmul(rand_mat.T, rand_mat)
        quad = sp.Quad(8, dim=dim, mean=mean, cov=cov)
        for i in range(len(cov)):
            mean_i = quad.integrate('v[{}]'.format(i))
            self.assertAlmostEqual(mean_i, mean[i])
            for j in range(len(cov)):
                fun = '(v[{}]-{})*(v[{}]-{})'.format(i, mean[i], j, mean[j])
                cov_ij = quad.integrate(fun)
                self.assertAlmostEqual(cov_ij, cov[i][j])

class TestHermiteTransform(unittest.TestCase):

    ## Test that hermite polynomials are orthonormal
    # def test_normalization_nodes(self):

