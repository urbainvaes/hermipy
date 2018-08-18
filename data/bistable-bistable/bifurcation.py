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

betas_up = np.load("white-noise-betas-up.npy")
betas_down = np.load("white-noise-betas-down.npy")
ms_up = np.load("white-noise-ms-up.npy")
ms_down = np.load("white-noise-ms-down.npy")
ax.plot(betas_up, ms_up, '0.6', linewidth=.5)
ax.plot(betas_down, ms_down, '0.6', linewidth=.5,
        label="White noise")

betas_up = np.load("epsilon=1o4-m0=.75-betas.npy")
betas_down = np.load("epsilon=1o4-m0=-.75-betas.npy")
ms_up = np.load("epsilon=1o4-m0=.75-ms.npy")
ms_down = np.load("epsilon=1o4-m0=-.75-ms.npy")
ax.plot(betas_up, ms_up, 'k.', markersize=.5)
ax.plot(betas_down, ms_down, 'k.', markersize=.5,
        label="$\\varepsilon = 1/4$")

betas_up = np.load("epsilon=1o10-m0=0.5-betas.npy")
betas_down = np.load("epsilon=1o10-m0=-0.5-betas.npy")
ms_up = np.load("epsilon=1o10-m0=0.5-ms.npy")
ms_down = np.load("epsilon=1o10-m0=-0.5-ms.npy")
ax.plot(betas_up, ms_up, 'b.', markersize=.5)
ax.plot(betas_down, ms_down, 'b.', markersize=.5,
        label="$\\varepsilon = 1/10$")

betas_up = np.load("epsilon=1o40-m0=0.5-betas.npy")
betas_down = np.load("epsilon=1o40-m0=-0.5-betas.npy")
ms_up = np.load("epsilon=1o40-m0=0.5-ms.npy")
ms_down = np.load("epsilon=1o40-m0=-0.5-ms.npy")
ax.plot(betas_up, ms_up, 'r.', markersize=.5)
ax.plot(betas_down, ms_down, 'r.', markersize=.5,
        label="$\\varepsilon = 1/40$")

ax.legend()
plt.savefig('full_bifurcation.eps', bbox_inches='tight')

plt.show()
