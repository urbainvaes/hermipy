import hermite.core as core
from scipy.special import binom
import numpy as np
import numpy.linalg as la


class Series:

    @staticmethod
    def natural_bissect(func, x1=0, x2=1000):
        f1, f2 = func(x1), func(x2)
        if f1 is 0:
            return x1
        elif f2 is 0:
            return x2
        assert f1*f2 < 0
        x3 = (x1+x2)//2
        f3 = func(x3)
        replace_arg = 'x2' if f1*f3 <= 0 else 'x1'
        new_args = {'x1': x1, 'x2': x2}
        new_args[replace_arg] = x3
        return Series.natural_bissect(func, **new_args)

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

        if degree is None:
            self.degree = Series.natural_bissect(
                    lambda x: int(binom(x + self.dim, x)) - len(self.coeffs))
        else:
            self.degree = degree

    def __add__(self, other):
        assert abs(self.dim - other.dim) < 1e-8
        assert la.norm(self.mean - other.mean, 2) < 1e-8
        assert la.norm(self.cov - other.cov, 2) < 1e-8
        new_coeffs = self.coeffs + other.coeffs
        return Series(new_coeffs, dim=self.dim, mean=self.mean, cov=self.cov)

    def project(self, direction):
        direction = core.to_numeric(direction)
        p_coeffs = core.project(self.coeffs, self.dim, direction)
        return Series(p_coeffs,
                      mean=[self.mean[direction]],
                      cov=[[self.cov[direction][direction]]])
