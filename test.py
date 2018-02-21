import numpy as np
import numpy.polynomial.hermite_e as herm
import numba as nb
import time

nodes_1d, weights_1d = herm.hermegauss(100)

@nb.jit(fastmath=True)
def function_to_integrate(x, y, z):
    return np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z) + np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z) + np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z) + np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z) + np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z)
    # return 1


@nb.njit(fastmath=True,parallel=True)
def integrate3(num_int_Points):
    result = 0.
    for i in range(num_int_Points):
        for j in range(num_int_Points):
            for k in range(num_int_Points):
                f_point = function_to_integrate(nodes_1d[i],nodes_1d[j],nodes_1d[k])
                weight_point = weights_1d[i] * weights_1d[j] * weights_1d[k]
                result += f_point * weight_point
    return result

start = time.time()

for i in range(100):
    result = integrate3(100)

print("Result: " + str(result))
print(time.time()-start)
