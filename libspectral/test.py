import spectral as sp
import numpy as np
import numba as nb
from libhermite import hermite as hm

function = 'v[0]'
nodes, weights = sp.hermegauss_nd([5, 3])

a=hm.double_cube()
a.append



# @nb.jit("f8(f8[:])", fastmath=True)
# def function(v):
#     return v[0]



# print(sp.integrate_with_quad(function, nodes, weights))

# print(sp.integrate(function, [100, 100, 100]))

# deg = [2**i for i in range(8)]
# for i in deg:
#     print(sp.integrate(lambda v: 1 + 0*v[0], [i]))
