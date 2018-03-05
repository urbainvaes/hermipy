import spectral as sp
import time
import numpy as np
from libhermite import hermite as hm

cube = hm.triple_products(5)

c=sp.convert_to_numpy_cube(hm.triple_products(3))

quad = sp.Quad(100)
quad.varf('1', 2)



# nodes, weights = sp.hermegauss_nd([5, 5, 5])
# function = 'cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2])'

# v = sp.transform_simple_quad('v[1]', 2, nodes, weights)
# print('python')
# for i in v: print(i)

# nodes, weights = sp.hermegauss_nd([100, 100, 100])
# function = 'cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2])'

# start = time.time()
# v = sp.transform_simple_quad(function, 20, nodes, weights)
# print("Time: " + str(time.time() - start))

# # print(v[0])

# start = time.time()
# for i in range(100):
#     result = sp.integrate_simple_quad(function, nodes, weights)
# print("Result: " + str(result) + ", Time: " + str(time.time() - start))

# # function = 'v[0]*v[1]'
# # start = time.time()
# # for i in range(100):
# #     result = sp.integrate_with_quad(function, nodes, weights)
# # print("Result: " + str(result) + ", Time: " + str(time.time() - start))
