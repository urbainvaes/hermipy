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

βmax, ε = 5.95, 0
betas = np.load("theta-1-white-noise-betas.npy")
ms = np.load("theta-1-white-noise-ms.npy")
condition = betas < βmax
betas, ms = np.extract(condition, betas), np.extract(condition, ms)
ax.plot(betas, ms, color=cmap(ε))
ax.plot(betas, -ms, color=cmap(ε), label="White noise")

ε = .1
betas = np.load("epsilon=1o10-betas.npy")
ms = np.load("epsilon=1o10-ms.npy")
condition = betas < βmax
betas, ms = np.extract(condition, betas), np.extract(condition, ms)
ax.plot(betas, ms, '.-', color=cmap(ε), markersize=.5)
ax.plot(betas, -ms, '.-', color=cmap(ε), markersize=.5,
        label="$\\varepsilon = "+str(ε)+"$")

ε = .2
betas = np.load("epsilon=1o5-betas.npy")
ms = np.load("epsilon=1o5-ms.npy")
condition = betas < βmax
betas, ms = np.extract(condition, betas), np.extract(condition, ms)
ax.plot(betas, ms, '.-', color=cmap(ε), markersize=.5)
ax.plot(betas, -ms, '.-', color=cmap(ε), markersize=.5,
        label="$\\varepsilon = "+str(ε)+"$")

# ε = .5
# betas = np.load("epsilon=1o2-m0=-.75-betas.npy")
# ms = np.load("epsilon=1o2-m0=-.75-ms.npy")
# condition = betas < βmax
# betas, ms = np.extract(condition, betas), np.extract(condition, ms)
# ax.plot([*betas, 0.6], [*ms, 0], '.-', color=cmap(ε), markersize=.5)
# ax.plot([*betas, 0.6], [*(-ms), 0], '.-', color=cmap(ε), markersize=.5,
#         label="$\\varepsilon = "+str(ε)+"$")

# βmax, ε = 2.5, 1
# betas = np.load("epsilon=1o1-m0=.75-betas.npy")
# ms = np.load("epsilon=1o1-m0=.75-ms.npy")
# condition = betas < βmax
# betas, ms = np.extract(condition, betas), np.extract(condition, ms)
# ax.plot(betas, ms, '.-', color=cmap(.999), markersize=.5)
# ax.plot(betas, -ms, '.-', color=cmap(.999), markersize=.5,
#         label="$\\varepsilon = "+str(ε)+"$")

ax.legend()
plt.savefig('full_bifurcation.eps', bbox_inches='tight')

plt.show()
