import hermite.core as core
import hermite.lib as lib
import hermite.position as pos

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
        mats = []
        for a in args:
            assert type(a) is Varf
            if type(a.matrix) is ss.csr_matrix:
                mats.append(a.matrix.todense().A)
            else:
                mats.append(a.matrix)
        tens_mat = core.tensorize(mats, sparse=sparse)
        tens_pos = pos.Position.tensorize([a.position for a in args])
        return Varf(tens_mat, tens_pos)

    def __init__(self, matrix, position, degree=None):
        self.matrix = matrix
        self.is_sparse = isinstance(matrix, ss.csr_matrix)
        self.position = position

        if degree is None:
            dim, npolys = self.position.dim, self.matrix.shape[0]
            self.degree = lib.natural_bissect(
                    lambda x: int(binom(x + dim, x)) - npolys)
        else:
            self.degree = degree

    def __eq__(self, other):
        assert type(other) is Varf
        return self.position == other.position \
            and la.norm(self.matrix - other.matrix) < very_small

    def __add__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_matrix = self.matrix + other

        elif type(other) is Varf:
            assert self.position == other.position
            new_matrix = self.matrix + other.matrix

        else:
            raise TypeError("Invalid type!)")

        return Varf(new_matrix, self.position)

    def __mul__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_matrix = self.matrix * other
            return Varf(new_matrix, self.position)

        elif type(other) is Varf:
            return Varf.tensorize([self, other])

        else:
            raise TypeError("Invalid type: " + str(type(other)))

    def project(self, direction):
        direction = core.to_numeric(direction)
        p_matrix = core.project(self.matrix, self.dim, direction)
        p_pos = self.position.project([direction])
        return Varf(p_matrix, p_pos)
