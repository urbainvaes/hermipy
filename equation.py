import sympy as sym
sym.init_printing()

# Space variables
x, y = sym.symbols('x y', real=True)

# Unknown function
f = sym.Function('f')(x, y)


def map_operator(operator, factor):
    return (operator.subs(f, f*factor).doit()/factor).simplify()


def forward_params():

    params, functions = {}, {}

    # Create missing real positive parameters
    options = {'real': True, 'positive': True}
    for param in ['βx', 'βy', 'γ', 'ε']:
        params[param] = sym.symbols(param, **options)

    # Create missing real parameters
    options = {'real': True}
    for param in ['θ', 'm']:
        params[param] = sym.symbols(param, **options)

    # Create missing functions
    for func, var in {'Vp': x, 'Vy': y}.items():
        functions[func] = sym.Function(func)(var)

    return {**params, **functions}


def forward(params):

    # Real parameters
    βx, βy, γ, ε, θ, m = (params[x] for x in ['βx', 'βy', 'γ', 'ε', 'θ', 'm'])

    # Functional parameters
    Vp, Vy = (params[x] for x in ['Vp', 'Vy'])

    d = sym.diff

    # Fokker planck operator
    return d(d(Vp, x)*f + θ*(x-m)*f + (1-γ)*sym.sqrt(2/βx)*y*f/ε, x) \
        + γ**2/βx * (1/ε) * d(d(f, x), x) \
        + (1/ε**2) * d(d(Vy, y) * f, y) \
        + (1/ε**2) * (1/βy) * d(d(f, y), y)


def solve_gaussian(operator):
    cx, cxy, cy = sym.symbols('sol_cx sol_cxy sol_cy')
    ansatz = sym.exp(-sym.Rational(1, 2)*(cx*x*x + 2*cxy*x*y + cy*y*y))
    res = (operator.subs(f, ansatz)/ansatz).doit().cancel()
    coeffs_dict = res.as_poly(x, y).as_dict()
    coeffs_list = list(coeffs_dict[k] for k in coeffs_dict)
    (sol_cx, sol_cxy, sol_cy) = sym.solve(coeffs_list, cx, cxy, cy)[0]
    solution = ansatz.subs([(cx, sol_cx), (cxy, sol_cxy), (cy, sol_cy)])
    assert operator.subs(f, solution).doit().expand() == 0
    determinant = sol_cx * sol_cy - sol_cxy * sol_cxy
    factor = sym.sqrt(determinant) / (2*sym.pi)
    solution = factor * solution
    return solution
