#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hermite_cpp as test
import numpy as np
import numpy.linalg as la
import time


def timeit(function):
    def wrapper(*args, **kwargs):
        tb = time.time()
        result = function(*args, **kwargs)
        te = time.time()
        print(te - tb)
        return result
    return wrapper

n = 10000

mat_cpp = timeit(test.test)(n)
mat_numpy = timeit(test.to_numpy)(mat_cpp)
# mat_cpp = timeit(test.to_mat)(mat_numpy)
bmat_cpp = timeit(test.to_bmat)(mat_numpy)
# mat_numpy1 = timeit(test.to_numpy)(mat_cpp)
mat_numpy2 = timeit(test.to_numpy)(bmat_cpp)

# print(np.max(mat_numpy1 - mat_numpy2))

# mat_cpp_contig = timeit(test.test3)(n)
# mat_numpy_contig = timeit(test.to_numpy)(mat_cpp_contig)
# mat_numpy_contig2 = timeit(test.test2)(n)

# print(la.norm(mat_numpy_contig - mat_numpy_contig2, 2))
# print(mat_numpy_contig)
# print(mat_numpy_contig2)


# print(mat
# print(mat_numpy)

# print(mat_usual.flags)
# print(mat_usual.strides)

# print(mat_usual)

# mat_contig = timeit(test.to_numpy)(mat_contig)
# print(mat_contig.flags)
# print(mat_contig.strides)

# print(mat_contig)

# print(la.norm(mat_usual - mat_contig, 2))
