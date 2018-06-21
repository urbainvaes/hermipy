#!/usr/bin/env python
# -*- coding: utf-8 -*-

import hermipy.quad as quad
import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2)

degree, n_points = 25, 200

q1 = quad.Quad.gauss_hermite(n_points, mean=[0], cov=[[0.1]])
s1 = q1.transform('cos(5*x)', degree) * q1.discretize(q1.weight())


q1.plot(s1, 1, ax1)
# s1.plot(ax1)

q2 = quad.Quad.gauss_hermite([n_points, n_points])
s2 = q2.transform('+ cos(x*x) + exp(x) + cos(x) * cos(y)', degree)
s2.plot(ax2)

plt.show()
