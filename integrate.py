import numpy as np
from libhermite import hermite

print(hermite.integrate(100, lambda x: np.sin(x)**2))
print(hermite.integrate_2d(100, lambda x, y: np.sin(x)**2))
print(hermite.integrate_3d(100, lambda x, y, z: np.sin(x)**2))
