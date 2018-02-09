import numpy as np
from libhermite import hermite

# print(hermite.integrate(64, lambda x: np.sin(x)**2))
# print(hermite.integrate_2d(64, lambda x, y: np.sin(x)**2))
# print(hermite.integrate_3d(64, lambda x, y, z: np.sin(x)**2))

print("Assemble quadratures...")
quad_64 = hermite.Quad(64, 3)
quad_smolyak = hermite.Quad(0, 3)

def function(x, y, z): return np.cos(2*x+2*y+2*z)
def poly(x, y, z): return x**2+y**3+(x+y+z)**4

print("Calculate 3D integral with Gauss-Hermite...")
result_64 = quad_64.integrate(function)
print("-> Result = " + str(result_64))

print("Calculate 3D integral with Smolyak...")
result_smolyak = quad_smolyak.integrate(function)
print("-> Result = " + str(result_smolyak))
