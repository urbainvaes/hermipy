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

    def inner(self, s2):
        assert type(self) is Series and type(s2) is Series
        assert self.degree == s2.degree
        d1, d2 = self.position.dirs, s2.position.dirs
        result = core.inner(self.coeffs, s2.coeffs, d1, d2)
        inner_pos = pos.Position.inner(self.position, s2.position)
        return Series(result, inner_pos, degree=self.degree)

    def project(self, directions):
        if type(directions) is int:
            directions = [directions]
        directions = core.to_numeric(directions)
        p_coeffs = core.project(self.coeffs, self.position.dim, directions)
        p_pos = self.position.project(directions)
        return Series(p_coeffs, p_pos)

    def subdegree(self, degree):
        assert degree <= self.degree
        n_polys = int(binom(degree + self.position.dim, degree))
        coeffs = self.coeffs[0:n_polys]
        return Series(coeffs, self.position, degree=degree)
