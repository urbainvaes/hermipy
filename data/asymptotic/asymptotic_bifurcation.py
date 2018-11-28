#!/usr/bin/env python

import numpy as np
import sympy as sym
import scipy.integrate as integrate
import scipy.optimize as optimize
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rc('font', size=18)
matplotlib.rc('font', family='serif')
matplotlib.rc('text', usetex=True)
x = sym.symbols('x')


εs_spectral_ou = [0.1, 0.2, 0.3, 0.4, 0.5, 1]
βs_spectral_ou = [2.152, 2.055, 1.908, 1.742, 1.57, 0.896]

εs_spectral_harmonic = [0.1, 0.2, 0.3, 0.4, 0.5]
βs_spectral_harmonic = [2.138, 2.127, 2.094, 2.012, 1.96]

# With correction for wrong effective drift
factor_degree_30 = 0.9728
βs_spectral_harmonic = np.asarray(βs_spectral_harmonic)/factor_degree_30


def Vt(m):
    θ = 1
    return x**4/4 - x**2/2 + θ*(x - m)**2/2


def moment(β, m, noise='ou'):

    rho_0 = sym.exp(-β*Vt(m))
    Z = integrate.quad(sym.lambdify(x, rho_0), -4, 4)[0] - 0
    rho_0 = rho_0 / Z

    C = sym.symbols('C')
    Vtm = Vt(m)

    # Power of first non-zero correction
    p = 2 if noise == 'ou' else 4

    if noise == 'ou':
        rho_p = C*β - 0.5*β*Vtm.diff(x)**2 + Vtm.diff(x, x)
        Cn_func = sym.solve(rho_p, C)[0]*rho_0
        Cn = integrate.quad(sym.lambdify(x, Cn_func), -4, 4)[0]
        rho_p = rho_0*rho_p.subs(C, Cn)

    elif noise == 'harmonic':
        dV, ddV, dddV = Vtm.diff(x), Vtm.diff(x, x), Vtm.diff(x, x, x)
        rho_p = (C*β - β*dV**2*ddV + β*sym.integrate(dV*ddV**2, x)
                 + 2*dV*dddV + ddV**2/2 - Vtm.diff(x, x, x, x)/β)
        Cn_func = sym.solve(rho_p, C)[0]*rho_0
        Cn = integrate.quad(sym.lambdify(x, Cn_func), -4, 4)[0]
        rho_p = rho_0*rho_p.subs(C, Cn)

    assert abs(integrate.quad(sym.lambdify(x, rho_0), -4, 4)[0] - 1) < 1e-8
    assert abs(integrate.quad(sym.lambdify(x, rho_p), -4, 4)[0] - 0) < 1e-8

    i0 = integrate.quad(sym.lambdify(x, x*rho_0), -4, 4)[0]
    ip = integrate.quad(sym.lambdify(x, x*rho_p), -4, 4)[0]

    return lambda ε: i0 + ε**p*ip


def critical(β, noise='ou'):
    εmin, εmax = 0, 1

    dm = .01
    moment_β = moment(β, dm, noise=noise)

    def objective(ε):
        return moment_β(ε)/dm - 1

    omin, omax = objective(εmin), objective(εmax)
    print('β = ' + str(β))

    if not omin < 0 or not omax > 0:
        return None

    return optimize.bisect(objective, εmin, εmax)


def bifurcation_data(noise):
    β, Δβ = .01, .03
    while critical(β, noise=noise) is None:
        β += Δβ

    εold, ε = - np.inf, critical(β, noise=noise)
    while ε > εold:
        βold, εold = β, ε
        β, ε = β + Δβ, critical(β + Δβ, noise=noise)

    Δε, εs, βs = 0.01, [], []
    while True:
        dεdβ = (ε - εold)/(β - βold)

        while True:
            Δβ = min(Δβ, .9*abs(Δε/dεdβ))
            newε = critical(β + Δβ, noise=noise)
            newε = newε if newε is not None else 0
            if abs(newε - ε) < Δε:
                break
            dεdβ = (newε - ε)/Δβ

        βold, β = β, β + Δβ
        εold, ε = ε, newε
        εs.append(ε), βs.append(β)

        if ε is 0:
            return βs, εs


try:
    βs_ou = np.load("betas-critical-ou.npy")
    βs_harmonic = np.load("betas-critical-harmonic.npy")
    εs_ou = np.load("epsilons-critical-ou.npy")
    εs_harmonic = np.load("epsilons-critical-harmonic.npy")
except IOError:
    βs_ou, εs_ou = bifurcation_data('ou')
    βs_harmonic, εs_harmonic = bifurcation_data('harmonic')

np.save("epsilons-critical-ou", εs_ou)
np.save("epsilons-critical-harmonic", εs_harmonic)
np.save("betas-critical-ou", βs_ou)
np.save("betas-critical-harmonic", βs_harmonic)

fig, ax = plt.subplots()
ax.plot(εs_ou, βs_ou, 'r-', label='Asymptotic - OU noise')
ax.plot(εs_harmonic, βs_harmonic, 'b-', label='Asymptotic - Harmonic noise')
ax.plot(εs_spectral_ou, βs_spectral_ou, 'k.', markersize=10,
        label='Spectral method - OU noise')
ax.plot(εs_spectral_harmonic, βs_spectral_harmonic, 'g.', markersize=10,
        label='Spectral method - Harmonic noise')
plt.legend()
ax.set_xlim((0, 1))
ax.set_xlabel('$\\varepsilon$')
ax.set_ylabel('$\\beta^c$')
plt.savefig("asymptotic_bifurction.eps", bbox_inches='tight')
plt.show()
