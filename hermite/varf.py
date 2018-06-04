import hermite.core as core
import hermite.lib as lib

from scipy.special import binom
import scipy.sparse as ss

import numpy as np
import numpy.linalg as la

import ipdb


very_small = 1e-10


class Varf:

    @staticmethod
    def tensorize(args, sparse=False):
        assert len(args) > 1
        dim, mean, cov, mats = 0, [], [], []
        for a in args:
            assert type(a) is Varf
            assert a.is_diag
            dim += a.dim
            mean.extend(a.mean)
            cov.extend(np.diag(a.cov))
            if type(a.matrix) is ss.csr_matrix:
                mats.append(a.matrix.todense().A)
            else:
                mats.append(a.matrix)
        mean = np.asarray(mean)
        cov = np.diag(cov)
        tens_mat = core.tensorize(mats, sparse=sparse)
        return Varf(tens_mat, dim=dim, mean=mean, cov=cov)

    def __init__(self, matrix, dim=1, mean=None, cov=None, degree=None):

        self.dim = dim
        self.mean = np.zeros(self.dim) if mean is None \
            else np.asarray(mean, float)
        self.cov = np.eye(self.dim) if cov is None \
            else np.asarray(cov, float)

        diag_cov = np.diag(np.diag(self.cov))
        self.is_diag = la.norm(self.cov - diag_cov, 2) < 1e-10

        self.matrix = matrix
        self.is_sparse = isinstance(matrix, ss.csr_matrix)

        if degree is None:
            npolys = self.matrix.shape[0]
            self.degree = lib.natural_bissect(
                    lambda x: int(binom(x + self.dim, x)) - npolys)
        else:
            self.degree = degree

    def __eq__(self, other):
        assert type(other) is Varf
        return self.dim == other.dim \
            and la.norm(self.mean - other.mean, 2) < very_small \
            and la.norm(self.cov - other.cov, 2) < very_small \
            and la.norm(self.matrix - other.matrix) < very_small

    def __add__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_matrix = self.matrix + other

        elif type(other) is Varf:
            assert self.dim == other.dim
            assert la.norm(self.mean - other.mean, 2) < very_small
            assert la.norm(self.cov - other.cov, 2) < very_small
            new_matrix = self.matrix + other.matrix

        else:
            raise TypeError("Invalid type!)")

        return Varf(new_matrix, dim=self.dim, mean=self.mean, cov=self.cov)

    def __mul__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_matrix = self.matrix * other
            return Varf(new_matrix, dim=self.dim, mean=self.mean, cov=self.cov)

        elif type(other) is Varf:
            return Varf.tensorize([self, other])

        else:
            raise TypeError("Invalid type: " + str(type(other)))

    def project(self, direction):
        direction = core.to_numeric(direction)
        p_matrix = core.project(self.matrix, self.dim, direction)
        return Varf(p_matrix,
                    mean=[self.mean[direction]],
                    cov=[[self.cov[direction][direction]]])
