import unittest
import numpy as np
import spectral as sp


class TestIntegrate(unittest.TestCase):

    def test_normalization_nodes(self):
        deg = [2**i for i in range(8)]
        for i in deg:
            integral = sp.integrate(lambda v: 1 + 0*v[0], [i])
            self.assertAlmostEqual(integral, 1)

    def test_normalization_dim(self):
        for i in range(1, 4):
            integral = sp.integrate(lambda v: 1 + 0*v[0], 100, dim=i)
            self.assertAlmostEqual(integral, 1)

    def test_mean(self):
        dim = 3
        mean = np.random.random(dim)
        for i in range(len(mean)):
            coord = sp.integrate(lambda v: v[i], 8, dim=dim, mean=mean)
            self.assertAlmostEqual(coord, mean[i])

    def test_covariance(self):
        dim = 3
        rand_mat = np.random.random((dim, dim))
        cov = np.matmul(rand_mat.T, rand_mat)
        for i in range(len(cov)):
            for j in range(len(cov)):
                cov_ij = sp.integrate(lambda v: v[i]*v[j], 8, dim=dim, cov=cov)
                self.assertAlmostEqual(cov_ij, cov[i][j])

    def test_quad(self):
        dim = 3
        rand_mat = np.random.random((dim, dim))
        mean = np.arange(dim)
        cov = np.matmul(rand_mat.T, rand_mat)
        quad = sp.Quad(8, dim=dim, mean=mean, cov=cov)
        for i in range(len(cov)):
            mean_i = quad.integrate(lambda v: v[i])
            self.assertAlmostEqual(mean_i, mean[i])
            for j in range(len(cov)):
                cov_ij = quad.integrate(lambda v: (v[i]-mean[i])*(v[j]-mean[j]))
                self.assertAlmostEqual(cov_ij, cov[i][j])
