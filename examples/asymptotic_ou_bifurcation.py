import sympy as sym
import hermipy as hm
import hermipy.equations as eq

sym.init_printing()
equation = eq.McKean_Vlasov
q = hm.Quad.gauss_hermite(200, dirs=[0], mean=[0], cov=[[.1]])
x = equation.x


def Vt(m):
    θ = 1
    return x**4/4 - x**2/2 + θ*(x - m)**2/2


def moment(β, Vt):

    ρ0 = sym.exp(-β*Vt)
    Z = q.integrate(ρ0, flat=True)
    ρ0 = ρ0 / Z

    C = sym.symbols('C')
    ρ2 = C*β - 0.5*β**Vt.diff(x)**2 + Vt.diff(x, x)
    Cn = q.integrate(sym.solve(ρ2, C)[0]*ρ0, flat=True)
    ρ2 = ρ0*ρ2.subs(C, Cn)

    assert abs(q.integrate(ρ0, flat=True) - 1) < 1e-8
    assert abs(q.integrate(ρ2, flat=True) - 0) < 1e-8

    return lambda ε: q.integrate(x*(ρ0 + ε**2*ρ2), flat=True)


def slope(β):
    dm = .01
    Vr, Vl = Vt(dm/2), Vt(-dm/2)
    mr, ml = moment(β, Vr), moment(β, Vl)
    return lambda ε: 1/dm*(mr(ε) - ml(ε))


def critical(ε):
    βmin, βmax = .5, 5

    def objective(β):
        return slope(β, ε) - 1

    assert objective(βmax) > 0
    assert objective(βmin) < 0

    β = (βmin + βmax)/2
    while True:
        fβ = objective(β)
        if abs(fβ) < 1e-8:
            return β

        if fβ > 0:
            βmax = β
        else:
            βmin = β

        print('β: ' + str(β) + ', f(β): ' + str(fβ))


critical(0.1)
