import numpy as np
import sympy as sy
import time

from libhermite import hermite

print("Assemble quadratures...")

my_quad = hermite.Quad(100, 3)
quad_smolyak = hermite.Quad(0, 3)

# Import attributes to python
weights = my_quad.weights
nodes = my_quad.nodes
np_weights=np.array(range(len(weights)))
np_nodes=np.ndarray((len(nodes),len(nodes[0])))
for i in range(len(weights)):
    np_weights[i] = weights[i]
    for j in range(len(nodes[0])):
        np_nodes[i][j] = nodes[i][j]


def function(x: float, y: float, z: float) -> float: return np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z) +np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z) +np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z) +np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z) +np.cos(2*x+2*y+2*z) + x*y + np.exp(-z*z)
function_string = 'cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2]) +cos(2*v[0]+2*v[1]+2*v[2]) + v[0]*v[1] + exp(-v[2]*v[2])'
def function(x, y, z): return x*x
function_string = 'v[0]*v[0]'
# def poly(x, y, z): return x**2+y**3+(x+y+z)**4

print("Calculate 3D integral with Gauss-Hermite...")
start = time.time()
result_cpp = my_quad.integrate(function)
end = time.time()
print("-> Result = " + str(result_cpp) + ", Time = " + str(end-start))

print("Calculate 3D integral with Gauss-Hermite using integrate_from_string...")
start = time.time()
result_cpp = my_quad.integrate_from_string(function_string)
end = time.time()
print("-> Result = " + str(result_cpp) + ", Time = " + str(end-start))

print("Calculate 3D integral with python loop...")
start = time.time()
result_python = 0.
for i in range(len(weights)):
    result_python += np_weights[i] * function(np_nodes[i][0], np_nodes[i][1], np_nodes[i][2])
end = time.time()
print("-> Result = " + str(result_python) + ", Time = " + str(end-start))


print("Calculate 3D integral with Smolyak...")
start = time.time()
result_smolyak = quad_smolyak.integrate(function)
end = time.time()
print("-> Result = " + str(result_smolyak) + ", Time = " + str(end-start))

