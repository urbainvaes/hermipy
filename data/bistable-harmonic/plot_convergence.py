import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rc('font', size=20)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)

degrees = np.load('convergence_degree_corrected_degrees.npy')
error_l1 = np.load('convergence_degree_corrected_error_l1.npy')
error_eig = np.load('convergence_degree_corrected_error_eig.npy')

degree = 50
fig, ax = plt.subplots()
xplot, yplot = degrees, error_l1
ax.semilogy(xplot, yplot, 'b.',
            label="$\\|\\rho_{{ {} }} - \\rho_d\\|_1$".format(degree))
coeffs = np.polyfit(xplot, np.log10(yplot), 1)
ax.semilogy(xplot, 10**coeffs[1] * 10**(coeffs[0]*xplot), 'b-')
xplot, yplot = degrees, error_eig
ax.semilogy(xplot, yplot, 'r.', label="$|\\lambda_0(d)|$")
coeffs = np.polyfit(xplot, np.log10(yplot), 1)
ax.semilogy(xplot, 10**coeffs[1] * 10**(coeffs[0]*xplot), 'r-')
ax.set_xlabel("$d$")
plt.legend(loc='upper right')
plt.savefig("errors.eps", bbox_inches='tight')
plt.show()
