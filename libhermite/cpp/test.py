import hermite_cpp as test
import numpy as np
import time


def timeit(function):
    def wrapper(*args, **kwargs):
        tb = time.time()
        result = function(*args, **kwargs)
        te = time.time()
        print(te - tb)
        return result
    return wrapper


A = timeit(test.test)(10000)
B = timeit(test.test3)(10000)
C = timeit(test.test2)(10000)

print(type(B))
A = timeit(test.to_numpy)(A)
B = timeit(test.to_numpy)(B)

print(A - C)
# print(type(B))
# print(B)

# print(B)
# C = timeit(np.array)(A)
# D = timeit(np.copy)(C)
