from .core import multi_indices

import sympy as sym
import math
import re


def stringify(function):
    if isinstance(function, int) or isinstance(function, float):
        return str(function)
    if isinstance(function, sym.Expr):
        function = sym.ccode(function)
    if isinstance(function, str):
        function = re.sub(r'\bx\b', 'v[0]', function)
        function = re.sub(r'\by\b', 'v[1]', function)
        function = re.sub(r'\bz\b', 'v[2]', function)
        function = re.sub(r'(?<=[v])([0-9]+)', r'[\1]', function)
    return function


def x_ify(function):
    if not isinstance(function, tuple(sym.core.all_classes)):
        return function
    symbols = list(function.free_symbols)
    if symbols == []:
        return function
    assert len(symbols) == 1
    return function.subs(symbols[0], sym.symbols('x'))


def split_product(expression, symbols):
    is_mul = isinstance(expression, sym.mul.Mul)
    args = expression.args if is_mul else [expression]
    result = {str(s): sym.Rational('1') for s in symbols}
    for arg in args:
        str_symbols = [str(s) for s in arg.free_symbols]
        if len(str_symbols) == 0:
            result['v[0]'] *= arg
        elif len(str_symbols) == 1:
            x_ified = stringify(str_symbols[0])
            result[x_ified] *= x_ify(arg)
        else:
            return False
    return result


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
