#!/usr/bin/env python

import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rc('font', size=14)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)

fig, ax = plt.subplots(1, 1)
ax.set_xlabel("$\\beta$")
ax.set_ylabel("$m$")

cmap = matplotlib.cm.get_cmap('viridis_r')
factor_degree_30 = 0.9728

βmin, βmax, ε = 1.6, 4.9, 0
betas = np.load("theta-1-white-noise-betas.npy")
ms = np.load("theta-1-white-noise-ms.npy")
condition = (betas > βmin) * (betas < βmax)
betas, ms = np.extract(condition, betas), np.extract(condition, ms)
ax.plot(betas, ms, color=cmap(ε))
ax.plot(betas, -ms, color=cmap(ε), label="White noise")

βmin, βmax, ε = 1.6, 5, .1
betas = np.load("deg=30-epsilon=1o10-betas.npy")/factor_degree_30
ms = np.load("deg=30-epsilon=1o10-ms.npy")
condition = (betas > βmin) * (betas < βmax)
betas, ms = np.extract(condition, betas), np.extract(condition, ms)
ax.plot(betas, ms, '.-', color=cmap(ε), markersize=.5)
ax.plot(betas, -ms, '.-', color=cmap(ε), markersize=.5,
        label="$\\varepsilon = "+str(ε)+"$")

βmin, βmax, ε = 1.6, 5, .2
betas = np.load("deg=30-epsilon=1o5-betas.npy")/factor_degree_30
ms = np.load("deg=30-epsilon=1o5-ms.npy")
condition = (betas > βmin) * (betas < βmax)
betas, ms = np.extract(condition, betas), np.extract(condition, ms)
ax.plot(betas, ms, '.-', color=cmap(ε), markersize=.5)
ax.plot(betas, -ms, '.-', color=cmap(ε), markersize=.5,
        label="$\\varepsilon = "+str(ε)+"$")


βmax, ε = 5, .5
betas = np.load("deg=30-epsilon=1o2-betas.npy")/factor_degree_30
ms = np.load("deg=30-epsilon=1o2-ms.npy")
condition = betas < βmax
betas, ms = np.extract(condition, betas), np.extract(condition, ms)
ax.plot(betas, ms, '.-', color=cmap(ε), markersize=.5)
ax.plot(betas, -ms, '.-', color=cmap(ε), markersize=.5,
        label="$\\varepsilon = "+str(ε)+"$")

ax.legend()
plt.savefig('full_bifurcation.eps', bbox_inches='tight')

plt.show()
