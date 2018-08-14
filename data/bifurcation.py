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

# remove_zoom = 48
# betas_zoom_up = np.load("betas-zoom-up.npy")[0:-remove_zoom]
# ms_zoom_up = np.load("ms-zoom-up.npy")[0:-remove_zoom]
# ax.plot(betas_zoom_up, ms_zoom_up, 'g.', markersize=.5)

# ε = 1/10
# betas_up = np.load("betas-up.npy")
# betas_down = np.load("betas-down.npy")
# ms_up = np.load("ms-up.npy")
# ms_down = np.load("ms-down.npy")
# ax.plot(betas_up, ms_up, 'g.', markersize=.5)
# ax.plot(betas_down, ms_down, 'g.', markersize=.5)

betas_up = np.load("white-noise-betas-up.npy")
betas_down = np.load("white-noise-betas-down.npy")
ms_up = np.load("white-noise-ms-up.npy")
ms_down = np.load("white-noise-ms-down.npy")
ax.plot(betas_up, ms_up, 'k-', linewidth=.5)
# ax.plot(betas_down, ms_down, 'k-', linewidth=.5)

betas_up = np.load("epsilon=1o20-betas-up.npy")
betas_down = np.load("epsilon=1o20-betas-down.npy")
ms_up = np.load("epsilon=1o20-ms-up.npy")
ms_down = np.load("epsilon=1o20-ms-down.npy")
ax.plot(betas_up, ms_up, 'b.', markersize=.5)
# ax.plot(betas_down, ms_down, 'b.', markersize=.5)

betas_up = np.load("epsilon=1o40-betas-up.npy")
betas_down = np.load("epsilon=1o40-betas-down.npy")
ms_up = np.load("epsilon=1o40-ms-up.npy")
ms_down = np.load("epsilon=1o40-ms-down.npy")
ax.plot(betas_up, ms_up, 'r.', markersize=.5)
# ax.plot(betas_down, ms_down, 'r.', markersize=.5)

plt.savefig('full_bifurcation.eps', bbox_inches='tight')

plt.show()
