from hermite.core import multi_indices

import sympy as sym
import math

def split_operator(op, func, order):
    variables = func.args
    result, rem, order = [], op.expand(), 2
    for m in multi_indices(len(variables), order):
        if rem == 0:
            result.append(0)
            continue
        test, der = 1, func
        for i, v in zip(m, variables):
            test *= v**i/math.factorial(i)
            der = sym.diff(der, v, i)
        remargs = rem.args if isinstance(rem, sym.add.Add) else [rem]
        term, rem = 0, 0
        for arg in remargs:  # Convoluted to avoid rounding errors
            termarg = arg.subs(func, test).doit()
            if termarg == 0:
                rem += arg
            else:
                term += termarg
        if isinstance(term, tuple(sym.core.all_classes)):
            term = sym.simplify(term)
        result.append(term)
    assert rem == 0
    return result
