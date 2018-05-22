# from hermite.core import multi_indices

import sympy as sym
from sympy.parsing.sympy_parser import parse_expr
import re

import pdb


class Function():

    # xyz notation
    xyz = list(sym.symbols('x y z v w', real=True))

    # subscript notation
    v_sub = [sym.symbols('v{}'.format(i), real=True)
             for i in range(len(xyz))]

    # array notation
    v_array = [sym.symbols('v[{}]'.format(i), real=True)
               for i in range(len(xyz))]

    # Conversion between them
    conv = {}
    for i in range(len(xyz)):
        conv[str(xyz[i])] = v_array[i]
        conv[str(v_sub[i])] = v_array[i]
        conv[str(v_array[i])] = v_array[i]

    def __init__(self, expr):

        if isinstance(expr, int) or isinstance(expr, float):
            expr = str(expr)

        if isinstance(expr, str):

            # Parser does not work with array notation
            for i in range(len(self.xyz)):
                expr = re.sub(r'\bv\[{}\]'.format(i), str(self.v_sub[i]), expr)

            expr = parse_expr(expr)

        for s in expr.free_symbols:

            # Ensure the resulting function is over reals
            if str(s) in self.conv:
                expr = expr.subs(s, self.conv[str(s)])

            else:
                raise ValueError("Unrecognized variable: " + str(s))

        self.sym_func = expr

    def __eq__(self, other):
        return self.sym_func == other.sym_func

    def __str__(self):
        return sym.ccode(self.sym_func)

    def as_string(self, format='array'):
        function = str(self)
        if format == 'array':
            return function
        if format == 'xyz':
            to = self.xyz
        elif format == 'sub':
            to = self.v_sub
        else:
            raise ValueError("Invalid format: " + format)
        for i in range(len(self.xyz)):
            function = re.sub(r'\bv\[{}\]'.format(i), str(to[i]), function)
        return function
