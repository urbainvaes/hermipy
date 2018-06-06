#  TODO: Add directions in Position (urbain, 06 Jun 2018)

import hermite.core as core
import hermite.lib as lib
import hermite.position as pos

from scipy.special import binom
import numpy as np
import numpy.linalg as la
import ipdb


very_small = 1e-10


class Series:

    @staticmethod
    def tensorize(args):
        assert len(args) > 1
        vecs = []
        for a in args:
            assert type(a) is Series
            vecs.append(a.coeffs)
        tens_vec = core.tensorize(vecs)
        tens_pos = pos.Position.tensorize([a.position for a in args])
        return Series(tens_vec, tens_pos)

    @staticmethod
    def inner(s1, s2):
        assert type(s1) is Series and type(s2) is Series
        assert s1.degree == s2.degree
        d1, d2 = s1.position.dirs, s2.position.dirs
        result = core.inner(s1.coeffs, s2.coeffs, d1, d2)
        inner_pos = pos.Position.inner(s1.position, s2.position)
        return Series(result, inner_pos, degree=s1)

    def __init__(self, coeffs, position, degree=None, norm=False):
        self.coeffs = coeffs/la.norm(coeffs, 2) if norm else coeffs
        self.position = position

        if degree is None:
            dim, npolys = self.position.dim, len(self.coeffs)
            self.degree = lib.natural_bissect(
                    lambda x: int(binom(x + dim, x)) - npolys)
        else:
            self.degree = degree

    def __eq__(self, other):
        assert type(other) is Series
        return self.position == other.position \
            and la.norm(self.coeffs - other.coeffs) < very_small

    def __add__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_coeffs = self.coeffs + other

        elif type(other) is Series:
            assert self.position == other.position
            new_coeffs = self.coeffs + other.coeffs

        else:
            raise TypeError("Invalid type!)")

        return Series(new_coeffs, self.position)

    def __mul__(self, other):

        if isinstance(other, (int, float, np.float64)):
            new_coeffs = self.coeffs * other
            return Series(new_coeffs, self.position)

        elif type(other) is Series:
            return Series.tensorize([self, other])

        else:
            raise TypeError("Invalid type: " + str(type(other)))

    def project(self, direction):
        direction = core.to_numeric(direction)
        p_coeffs = core.project(self.coeffs, self.position.dim, direction)
        p_pos = self.position.project([direction])
        return Series(p_coeffs, p_pos)
