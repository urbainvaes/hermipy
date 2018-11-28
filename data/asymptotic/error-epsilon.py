import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rc('font', size=20)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)

epsilon, v, errors = [], [], []
for i in range(0, 6):
    epsilon.append(1/2**i)
    v.append(np.load("errors-epsilon-{}.npy".format(epsilon[-1])))
    errors.append(v[-1][-1])
epsilon, errors = np.array(epsilon), np.array(errors)

for i, y in enumerate(v):
    plt.semilogy(range(5, 5 + len(y)), y, '.', basey=2, 
                 label='$\\varepsilon = {}$'.format(epsilon[i]))
plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.xlabel("Degree of approximation")
plt.savefig('degree-difference.eps', bbox_inches='tight')
plt.show()

coeffs = np.polyfit(np.log2(epsilon), np.log2(errors), 1)
plt.loglog(epsilon, 2**coeffs[1] * epsilon**coeffs[0], 'r-',
            label='$y = {:.2f} \\times 2^{{ {:.2f} \\, x}}$'.
                  format(2**coeffs[1], coeffs[0]), basey=2)
plt.loglog(epsilon, errors, '.', basey=2, basex=2, markersize=10)
plt.legend()
plt.xlabel("$\\varepsilon$")
plt.savefig('epsilon-error.eps', bbox_inches='tight')
plt.show()
