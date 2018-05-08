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

n = 5

mat_cpp = timeit(test.test)(n)
mat_numpy = timeit(test.to_numpy)(mat_cpp)
mat_cpp = timeit(test.to_cpp)(mat_numpy)
mat_numpy = timeit(test.to_numpy)(mat_cpp)
mat_cpp = timeit(test.to_cpp)(mat_numpy.T)

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
