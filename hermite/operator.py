from hermite.core import multi_indices
import sympy as sym
import math


class Operator():

    def __init__(self, operator, function):
        self.function = function
        self.variables = function.args
        self.dimension = len(self.variables)
        self.operator = operator.expand()

    def split(self, order):
        result, rem, order = [], self.operator, 2
        for m in multi_indices(self.dimension, order):
            if rem == 0:
                result.append(0)
                continue
            test, der = 1, self.function
            for i, v in zip(m, self.variables):
                test *= v**i/math.factorial(i)
                der = sym.diff(der, v, i)
            remargs = rem.args if isinstance(rem, sym.add.Add) else [rem]
            term, rem = 0, 0
            for arg in remargs:  # Convoluted to avoid rounding errors
                termarg = arg.subs(self.function, test).doit()
                if termarg == 0:
                    rem += arg
                else:
                    term += termarg
            if isinstance(term, tuple(sym.core.all_classes)):
                term = sym.simplify(term)
            result.append(term)
        assert rem == 0
        return result
