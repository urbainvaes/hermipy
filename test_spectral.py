import unittest
import numpy as np
import spectral as sp
import numpy.polynomial.hermite_e as herm
import math

import importlib
importlib.reload(sp)


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

    def test_constant(self):
        degree = 30
        n_points = [degree, degree, degree]
        nodes, weights = sp.hermegauss_nd(n_points)
        coeffs = sp.transform_simple_quad('1', degree, nodes, weights)
        for i in range(len(coeffs)):
            target_value = 1. if i == 0 else 0.
            self.assertAlmostEqual(coeffs[i], target_value)

    def test_consistent_eval(self):
        degree = 10
        n_points = 10
        nodes, weights = sp.hermegauss_nd(n_points)
        nodes_scipy, weights_scipy = herm.hermegauss(n_points)
        for i in range(degree):
            coeffs = np.zeros(degree + 1)
            coeffs[i] = 1
            factor = math.sqrt(math.factorial(i))
            hi_nodes_scipy = herm.hermeval(nodes_scipy, coeffs)
            hi_nodes = sp.eval_simple_quad(coeffs, degree, nodes)
            diff = sum(abs(hi_nodes_scipy/factor - hi_nodes))
            self.assertAlmostEqual(diff, 0)

    def test_forward_backward(self):
        degree = 10
        nodes, weights = sp.hermegauss_nd(degree + 1)
        f_hermite = np.random.random(degree + 1)
        f_grid = sp.eval_simple_quad(f_hermite, degree, nodes)
        f_hermite_new = sp.transform_simple_quad(f_grid, degree, nodes, weights)
        diff = sum(abs(f_hermite - f_hermite_new))
        self.assertAlmostEqual(diff, 0)

    # def test_quad_transform(self):
        # dim = 3
        # n_points 
        # rand_mat = np.random.random((dim, dim))
        # mean = np.random.random(dim)
        # cov = np.matmul(rand_mat.T, rand_mat)
        # quad = sp.Quad(8, dim=dim, mean=mean, cov=cov)
        # n_points = [degree, degree, degree]
        # nodes, weights = sp.hermegauss_nd(n_points)
        # coeffs = sp.transform_simple_quad('1', degree, nodes, weights)
        # for i in range(len(coeffs)):
        #     target_value = 1. if i == 0 else 0.
        #     self.assertAlmostEqual(coeffs[i], target_value)

