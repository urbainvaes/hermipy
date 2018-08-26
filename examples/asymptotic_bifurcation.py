import multiprocessing
import numpy as np
import sympy as sym
import hermipy as hm
import hermipy.equations as eq
import scipy.integrate as integrate
import matplotlib
import matplotlib.pyplot as plt

sym.init_printing()
equation = eq.McKean_Vlasov
q = hm.Quad.gauss_hermite(200, dirs=[0], mean=[0], cov=[[.1]])
x = equation.x


def Vt(m):
    θ = 1
    return x**4/4 - x**2/2 + θ*(x - m)**2/2


def moment(β, m, noise='ou'):

    rho_0 = sym.exp(-β*Vt(m))
    Z = q.integrate(rho_0, flat=True)
    Z = integrate.quad(sym.lambdify(x, rho_0), -4, 4)[0] - 0
    rho_0 = rho_0 / Z

    C = sym.symbols('C')
    Vtm = Vt(m)

    # Power of first non-zero correction
    p = 2 if noise == 'ou' else 4

    if noise == 'ou':
        rho_p = C*β - 0.5*β*Vtm.diff(x)**2 + Vtm.diff(x, x)
        Cn = q.integrate(sym.solve(rho_p, C)[0]*rho_0, flat=True)
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

    # Check with hermipy integration
    check_consistency_hermipy = False
    if check_consistency_hermipy:
        hermipy_err_1 = abs(q.integrate(rho_0, flat=True) - 1)
        hermipy_err_2 = abs(q.integrate(rho_p, flat=True) - 0)
        if hermipy_err_1 > 1e-8 or hermipy_err_2 > 1e-8:
            print("Warning: hermipy gives different result: ",
                  + hermipy_err_1, hermipy_err_2)

    i0 = integrate.quad(sym.lambdify(x, x*rho_0), -4, 4)[0]
    ip = integrate.quad(sym.lambdify(x, x*rho_p), -4, 4)[0]

    return lambda ε: q.integrate(i0 + ε**p*ip)


def slope(β, noise='ou'):
    dm = .01
    return lambda ε: moment(β, dm, noise=noise)(ε)/dm


def critical(β, noise='ou'):
    εmin, εmax = 0, 1

    dm = .01
    moment_β = moment(β, dm, noise=noise)

    def objective(ε):
        return moment_β(ε)/dm - 1

    if not objective(εmin) < 0 and objective(εmax) > 0:
        return None

    while True:
        ε = (εmin + εmax)/2
        fε = objective(ε)

        if abs(fε) < 1e-8:
            return ε

        if fε > 0:
            εmax = ε
        else:
            εmin = ε

        print('β: ' + str(β) + ', ε: ' + str(ε) + ', f(ε): ' + str(fε))


βs = np.arange(.2, 2.5, .02)
εs_ou = [critical(β, noise='ou') for β in βs]
εs_harmonic = [critical(β, noise='harmonic') for β in βs]

condition = [ε is not None for ε in εs_ou]
βs_ou = np.extract(condition, βs)
εs_ou = np.extract(condition, εs_ou)

condition = [ε is not None for ε in εs_harmonic]
βs_harmonic = np.extract(condition, βs)
εs_harmonic = np.extract(condition, εs_harmonic)

np.save("epsilons-critical-ou", εs_ou)
np.save("epsilons-critical-harmonic", εs_harmonic)
np.save("betas-critical-ou", βs_ou)
np.save("betas-critical-harmonic", βs_harmonic)

fig, ax = plt.subplots()
ax.plot(εs_ou, βs_ou, 'r.', label='OU noise')
ax.plot(εs_harmonic, βs_harmonic, 'b+', label='Harmonic noise')
plt.legend()
plt.show()

import ipdb; ipdb.set_trace()
