# Copyright (C) 2018 Urbain Vaes
#
# This file is part of hermipy, a python/C++ library for automating the
# Hermite Galerkin method.
#
# hermipy is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# hermipy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import sympy as sym
from sympy.parsing.sympy_parser import parse_expr
import re


class Function():

    # xyz notation
    xyz = list(sym.symbols('x y z v w', real=True))

    # subscript notation
    v_sub = [sym.symbols('v{}'.format(i), real=True)
             for i in range(len(xyz))]

    # array notation
    v_array = [sym.symbols('v[{}]'.format(i), real=True)
               for i in range(len(xyz))]

    conv = {}
    for i in range(len(xyz)):
        conv[str(xyz[i])] = v_array[i]
        conv[str(v_sub[i])] = v_array[i]
        conv[str(v_array[i])] = v_array[i]

    def __init__(self, expr, dim=None, dirs=None, allow_sym=False):

        self.allow_sym = allow_sym

        if isinstance(expr, int) or isinstance(expr, float):
            expr = str(expr)

        if isinstance(expr, str):

            # Parser does not work with array notation
            for i in range(len(self.xyz)):
                expr = re.sub(r'\bv\[{}\]'.format(i), str(self.v_sub[i]), expr)

            expr = parse_expr(expr)

        assert isinstance(expr, tuple(sym.core.all_classes))
        for s in expr.free_symbols:

            # Ensure the resulting function is over reals
            if str(s) in self.conv:
                expr = expr.subs(s, self.conv[str(s)])

            elif not allow_sym:
                raise ValueError("Unrecognized variable: " + str(s))

        self.sym_func = expr

        if dirs is not None:
            self.dirs, self.dim = dirs, len(dirs)
            assert self.dirs == sorted(self.dirs)
            assert dim is None or dim == self.dim
        else:
            if dim is not None:
                self.dim = dim
            else:
                for i in range(len(self.v_array)):
                    if self.v_array[- 1 - i] in expr.free_symbols:
                        self.dim = len(self.v_array) - i
                        break
                else:
                    self.dim = 0
            self.dirs = list(range(self.dim))

        if not allow_sym:
            if len(expr.free_symbols) > self.dim:
                raise ValueError("Too many free symbols: " + str(expr)
                                 + ", dim: " + str(self.dim))

            symbols = [self.v_array[d] for d in self.dirs]
            for s in expr.free_symbols:
                if s not in symbols:
                    print("Unknown symbol!")
                    import pdb
                    pdb.set_trace()

    def __eq__(self, other):
        return self.dirs == other.dirs and \
            self.sym_func == other.sym_func

    def __str__(self):
        return str(self.sym_func)

    def __repr__(self):
        variables = [self.v_array[d] for d in self.dirs]
        return str(variables) + " --> " + str(self)

    def __mul__(self, other):
        if type(other) in (float, int):
            new_sym_func = self.sym_func*other
        elif type(other) is Function:
            assert self.dirs == other.dirs
            new_sym_func = self.sym_func * other.sym_func
        else:
            raise ValueError("Unsupported type: " + str(type(other)))
        return Function(new_sym_func, dirs=self.dirs)

    def __add__(self, other):
        assert self.dirs == other.dirs
        return Function(self.sym_func + other.sym_func,
                        dirs=self.dirs)

    def as_format(self, format):
        if format == 'array':
            return self.sym_func

        if format == 'xyz':
            to = self.xyz
        elif format == 'sub':
            to = self.v_sub
        else:
            raise ValueError("Invalid format: " + format)

        result = self.sym_func
        for i in range(len(self.xyz)):
            result = result.subs(self.v_array[i], to[i])

        return result

    def as_string(self, format='array', toC=False):
        if toC:
            # Convert dirs[i] -> i, to ease discretization
            function = sym.ccode(self.sym_func)
            for i in range(self.dim):
                if self.dirs[i] is not i:
                    assert not re.search(r'\bv\[{}\]'.format(i), function)
                    function = re.sub(r'\bv\[{}\]'.format(self.dirs[i]),
                                      'v[{}]'.format(i), function)
            return function
        return sym.ccode(self.as_format(format))

    def with_vars(self, variables):
        sym_func = self.sym_func
        for i, d in enumerate(self.dirs):
            sym_func = sym_func.subs(self.v_array[d], variables[i])
        return sym_func

    def split(self):
        variables = [self.v_array[d] for d in self.dirs]
        is_add = isinstance(self.sym_func, sym.add.Add)
        add_terms = self.sym_func.args if is_add else [self.sym_func]
        to_return = []
        for term in add_terms:
            is_mul = isinstance(term, sym.mul.Mul)
            mul_terms = term.args if is_mul else [term]
            result, tensorizable = [sym.Rational('1')] * (self.dim + 1), True
            for arg in mul_terms:
                vars_arg = [v for v in variables if v in arg.free_symbols]
                if len(vars_arg) is 0:
                    result[-1] *= (arg)  # Store constant in last element
                else:
                    if len(vars_arg) is 1:
                        symbol = vars_arg[0]
                        result[variables.index(symbol)] *= arg
                    else:
                        result[0] *= arg
                        tensorizable = False
            if tensorizable:
                for i in range(self.dim):
                    result[i] = Function(result[i], dirs=[self.dirs[i]],
                                         allow_sym=self.allow_sym)
            else:
                func = Function(term/result[-1], dirs=self.dirs,
                                allow_sym=self.allow_sym)
                result = [func, result[-1]]
            to_return.append(result)

        # factorized = []
        # for i, term in enumerate(to_return):
        #     for j in range(i, len(to_return)):
        #         other = to_return[j]
        #         if len(term

        return to_return

    @staticmethod
    def sanitize(expr, max_denom=1e8):
        result, expr = sym.Integer(0), (sym.Integer(0) + expr).expand()
        is_add = isinstance(expr, sym.add.Add)
        add_terms = expr.args if is_add else [expr]
        for aterm in add_terms:
            result_tmp = sym.Integer(1)
            is_mul = isinstance(aterm, sym.mul.Mul)
            mul_terms = aterm.args if is_mul else [aterm]
            for mterm in mul_terms:
                if type(mterm) is sym.numbers.Float:
                    if abs(mterm) < 1e-14:
                        mterm = 0
                    mterm = sym.Rational(mterm).limit_denominator(max_denom)
                result_tmp *= mterm
            result += result_tmp
        return result
