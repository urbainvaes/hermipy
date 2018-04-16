import sympy as sym
import sympy.printing as syp
sym.init_printing()

# Space variables
x, y = sym.symbols('x y', real=True)

# Unknown function
f =  sym.Function('f')(x, y)

def fill_params(params, functions):

    # Create missing real positive parameters
    options = { 'real' : True, 'positive' : True }
    for param in ['βx', 'βy', 'γ', 'ε']:
        if param not in params:
            params[param] = sym.symbols(param, **options)

    # Create missing real parameters
    options = { 'real' : True }
    for param in ['θ', 'm']:
        if param not in params:
            params[param] = sym.symbols(param, **options)

    # Create missing functions
    for fx in ['Vp', 'Vq']:
        if fx not in functions:
            functions[fx] = sym.Function(fx)(x)

    if 'Vy' not in functions:
        functions['Vy'] = sym.Function('Vy')(y)


def operators(params, functions, λ = sym.Rational(.5)):

    # Parameter of the mapping
    λ = sym.Rational(1,2)

    # Real parameters
    βx, βy, γ, ε, θ, m = (params[x] for x in ['βx', 'βy', 'γ', 'ε', 'θ', 'm'])

    # Functional parameters
    Vp, Vq, Vy = (functions[x] for x in ['Vp', 'Vq', 'Vy'])

    # Factor in mapping
    factor_x = sym.exp(-βx*(λ*Vq + (1-λ)*Vp))
    factor_y = sym.exp(-βy*Vy)
    factor = factor_x * factor_y

    d = sym.diff

    # Fokker planck operator
    forward = d(d(Vp, x)*f + θ*(x-m)*f + (1-γ)*sym.sqrt(2/βx)*y*f/ε, x) \
              + γ * (1/ε) * d(d(f, x), x) \
              + (1/ε**2) * d(d(Vy, y) * f, y) \
              + (1/ε**2) * (1/βy) * d(d(f, y), y)

    # Mapped operator, with solution in Gaussian-weighted space
    backward = (forward.subs(f, f*factor).doit()/factor).simplify()

    return forward, backward, factor_x, factor_y, factor


def solve_gaussian(operator):
    cx, cxy, cy = sym.symbols('cx cxy cy')
    ansatz = sym.exp(-(cx*x*x + cxy*x*y + cy*y*y))
    res = (operator.subs(f,ansatz)/ansatz).doit().cancel()
    coeffs_dict = res.as_poly(x, y).as_dict()
    coeffs_list = list(coeffs_dict[k] for k in coeffs_dict)
    (sol_cx, sol_cxy, sol_cy) = sym.solve(coeffs_list, cx, cxy, cy)[0]
    solution = ansatz.subs([(cx, sol_cx), (cxy, sol_cxy), (cy, sol_cy)])
    assert operator.subs(f,solution).doit().expand() == 0
    return solution
