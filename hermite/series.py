import hermite.core as core
import hermite.lib as lib

from scipy.special import binom
import numpy as np
import numpy.linalg as la
import ipdb


very_small = 1e-10


class Series:

    @staticmethod
    def tensorize(args):
        assert len(args) > 1
        dim, mean, cov, vecs = 0, [], [], []
        for a in args:
            assert type(a) is Series
            assert a.is_diag
            dim += a.dim
            mean.extend(a.mean)
            cov.extend(np.diag(a.cov))
            vecs.append(a.coeffs)
        mean = np.asarray(mean)
        cov = np.diag(cov)
        tens_vec = core.tensorize(vecs)
        return Series(tens_vec, dim=dim, mean=mean, cov=cov)

    def __init__(self, coeffs, dim=1, mean=None, cov=None,
                 degree=None, norm=False):
        self.coeffs = coeffs/la.norm(coeffs, 2) if norm else coeffs

        self.dim = dim
        self.mean = np.zeros(self.dim) if mean is None \
            else np.asarray(mean, float)
        self.cov = np.eye(self.dim) if cov is None \
            else np.asarray(cov, float)

        eigval, eigvec = la.eig(self.cov)
        self.factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

        diag_cov = np.diag(np.diag(self.cov))
        self.is_diag = la.norm(self.cov - diag_cov, 2) < 1e-10

        if degree is None:
            self.degree = lib.natural_bissect(
                    lambda x: int(binom(x + self.dim, x)) - len(self.coeffs))
        else:
            self.degree = degree

    def __eq__(self, other):
        assert type(other) is Series
        return self.dim == other.dim \
            and la.norm(self.mean - other.mean, 2) < very_small \
            and la.norm(self.cov - other.cov, 2) < very_small \
            and la.norm(self.coeffs - other.coeffs) < very_small

    def __add__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_coeffs = self.coeffs + other

        elif type(other) is Series:
            assert self.dim == other.dim
            assert la.norm(self.mean - other.mean, 2) < very_small
            assert la.norm(self.cov - other.cov, 2) < very_small
            new_coeffs = self.coeffs + other.coeffs

        else:
            raise TypeError("Invalid type!)")

        return Series(new_coeffs, dim=self.dim, mean=self.mean, cov=self.cov)

    def __mul__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_coeffs = self.coeffs * other
            return Series(new_coeffs, dim=self.dim,
                          mean=self.mean, cov=self.cov)

        elif type(other) is Series:
            return Series.tensorize([self, other])

        else:
            raise TypeError("Invalid type: " + str(type(other)))

    def project(self, direction):
        direction = core.to_numeric(direction)
        p_coeffs = core.project(self.coeffs, self.dim, direction)
        return Series(p_coeffs,
                      mean=[self.mean[direction]],
                      cov=[[self.cov[direction][direction]]])
