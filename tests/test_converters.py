import hermite.core as core
import hermite_cpp as hm
import unittest

import numpy as np
import numpy.linalg as la
import scipy.sparse as ss


class TestConversions(unittest.TestCase):

    def test_simple_sparse(self):
        matrix = [[0, 1, 0, 0, 0],
                  [0, 2, 0, 3, 0],
                  [0, 0, 4, 0, 0],
                  [5, 6, 0, 0, 0],
                  [7, 0, 0, 0, 8]]
        matrix = np.array(matrix, dtype=float)
        csr_mat = ss.csr_matrix(matrix)
        spmat = core.convert_to_cpp_sparse(csr_mat)
        cpp_mat = hm.full(spmat)
        self.assertTrue(la.norm(hm.to_numpy(cpp_mat) - matrix) < 1e-10)
