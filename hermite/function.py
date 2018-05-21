# from hermite.core import multi_indices

import sympy as sym
from sympy.parsing.sympy_parser import parse_expr

import re


class Function():

    xyz = list(sym.symbols('x y z w'))
    v = [sym.symbols('v[{}]'.format(i)) for i in range(len(xyz))]

    conv_xyz = {xyz[i]: v[i] for i in range(len(xyz))}

    def __init__(self, expr):

        # Support for v0 v1 v2 notation
        expr = re.sub(r'(?<=[v])([0-9]+)', r'[\1]', expr)

        if isinstance(expr, int) or isinstance(expr, float):
            expr = str(expr)
        if isinstance(expr, str):
            sym_func = parse_expr(expr)
        for s in sym_func.free_symbols:
            if s in self.conv_xyz:
                sym_func = sym_func.subs(s, self.conv_xyz[s])
            else:
                raise ValueError("Unrecognized variable")
        self.sym_func = sym_func

    def __eq__(self, other):
        return self.sym_func == other.sym_func

    def __str__(self):
        return sym.ccode(self.sym_func)

    def print(self, format='array'):
        function = self.__str__()
        if format == 'array':
            return function
        elif format == 'xyz':
            function = re.sub(r'\bv[0]', 'x', function)
            function = re.sub(r'\bv[1]', 'y', function)
            function = re.sub(r'\bv[2]', 'z', function)
        return function
