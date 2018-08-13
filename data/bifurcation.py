import numpy as np
import matplotlib
matplotlib.use('Agg')
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

# remove = 30
# betas_up = np.load("betas-up.npy")[0:-remove]
# betas_down = np.load("betas-down.npy")[0:-remove]
# ms_up = np.load("ms-up.npy")[0:-remove]
# ms_down = np.load("ms-down.npy")[0:-remove]
# ax.plot(betas_up, ms_up, 'r.', markersize=.5)
# ax.plot(betas_down, ms_down, 'b.', markersize=.5)

betas_up_small_epsilon = np.load("betas-other-epsilon-up.npy")
betas_down_small_epsilon = np.load("betas-other-epsilon-down.npy")
ms_up_small_epsilon = np.load("ms-other-epsilon-up.npy")
ms_down_small_epsilon = np.load("ms-other-epsilon-down.npy")
ax.plot(betas_up_small_epsilon, ms_up_small_epsilon, 'k.', markersize=.5)
ax.plot(betas_down_small_epsilon, ms_down_small_epsilon, 'k.', markersize=.5)

plt.savefig('full_bifurcation.eps', bbox_inches='tight')

plt.show()
