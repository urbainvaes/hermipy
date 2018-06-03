import numpy.polynomial.hermite_e as herm
import numpy as np

import sympy as sym
import math

from hermite.core import multi_indices


def hermegauss_nd(n_points):
    dim = len(n_points)
    nodes_multidim = []
    weights_multidim = []
    for i in range(dim):
        nodes_1d, weights_1d = herm.hermegauss(n_points[i])
        weights_1d = weights_1d/np.sqrt(2*np.pi)  # Normalize
        nodes_multidim.append(nodes_1d)
        weights_multidim.append(weights_1d)
    return nodes_multidim, weights_multidim


def natural_bissect(func, x1=0, x2=1000):
    f1, f2 = func(x1), func(x2)
    if f1 is 0:
        return x1
    elif f2 is 0:
        return x2
    assert f1*f2 < 0
    x3 = (x1+x2)//2
    f3 = func(x3)
    replace_arg = 'x2' if f1*f3 <= 0 else 'x1'
    new_args = {'x1': x1, 'x2': x2}
    new_args[replace_arg] = x3
    return natural_bissect(func, **new_args)


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
