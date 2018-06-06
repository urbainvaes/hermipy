import numpy as np
import numpy.linalg as la
import sympy as sym


class Position:

    very_small = 1e-10

    @staticmethod
    def tensorize(args):
        dim, mean, cov = 0, [], []
        for a in args:
            assert type(a) is Position
            assert a.is_diag
            dim += a.dim
            mean.extend(a.mean)
            cov.extend(np.diag(a.cov))
        return Position(dim=dim, mean=mean, cov=np.diag(cov))

    def __init__(self, dim=None, mean=None, cov=None, dirs=None):

        if mean is not None:
            self.dim = len(mean)
        elif cov is not None:
            self.dim = len(cov)
        else:
            assert dim is not None
            self.dim = dim

        self.mean = np.zeros(self.dim) if mean is None \
            else np.asarray(mean, float)

        self.cov = np.eye(self.dim) if cov is None \
            else np.asarray(cov, float)

        self.dirs = np.arange(dim) if dirs is None \
            else np.asarray(dirs)

        eigval, eigvec = la.eig(self.cov)
        self.factor = np.matmul(eigvec, np.sqrt(np.diag(eigval)))

        diag_cov = np.diag(np.diag(self.cov))
        self.is_diag = la.norm(self.cov - diag_cov, 2) < 1e-10

    def __eq__(self, other):
        return self.dim == other.dim \
            and la.norm(self.mean - other.mean, 2) < self.very_small \
            and la.norm(self.cov - other.cov, 2) < self.very_small

    def __mul__(self, other):
        assert type(other) is Position
        return Position.tensorize([self, other])

    def __hash__(self):
        return hash(frozenset({
            self.dim,
            hash(frozenset(self.mean.flatten())),
            hash(frozenset(self.cov.flatten()))}))

    def weight(self):
        var = [sym.symbols('v' + str(i), real=True) for i in range(self.dim)]
        inv_cov = la.inv(self.cov)
        potential = 0.5 * inv_cov.dot(var - self.mean).dot(var - self.mean)
        normalization = 1/(np.sqrt((2*np.pi)**self.dim * la.det(self.cov)))
        return normalization * sym.exp(-potential)

    def project(self, directions):
        assert self.is_diag
        dim = len(directions)
        mean = np.zeros(dim)
        cov = np.zeros((dim, dim))
        for i in range(dim):
            d = directions[i]
            assert d < self.dim
            mean[i] = self.mean[d]
            cov[i][i] = self.cov[d][d]
        return Position(dim=dim, mean=mean, cov=cov)
