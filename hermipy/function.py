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
import re


class Function():

    # xyz notation
    xyz = list(sym.symbols('x y z v w', real=True))

    # subscript notation
    x_sub = [sym.symbols('x{}'.format(i), real=True)
             for i in range(len(xyz))]

    conv = {}
    for i in range(len(xyz)):
        conv[str(xyz[i])] = x_sub[i]
        conv[str(x_sub[i])] = x_sub[i]

    def __init__(self, expr, dim=None, dirs=None, allow_sym=False):

        self.allow_sym = allow_sym

        # So that symbols are recognized despite possibly different properties
        expr = str(expr)

        if isinstance(expr, str):
            expr = sym.sympify(expr, locals=self.conv)

        assert isinstance(expr, tuple(sym.core.all_classes))
        for s in expr.free_symbols:
            if s not in self.x_sub and not allow_sym:
                raise ValueError("Unrecognized variable: " + str(s))

        self.sym_func = expr

        variables = expr.free_symbols.intersection(self.x_sub)
        if dirs is not None:
            self.dirs, self.dim = dirs, len(dirs)
            assert self.dirs == sorted(self.dirs)
            assert dim is None or dim == self.dim
            assert variables.issubset({self.x_sub[d] for d in self.dirs})
        else:
            # We do not minimize the number of dimensions
            self.dim = max(self.x_sub.index(s) for s in variables) + 1
            self.dirs = list(range(self.dim))

    def __eq__(self, other):
        return self.dirs == other.dirs and \
            self.sym_func == other.sym_func

    def __repr__(self):
        variables = [self.x_sub[d] for d in self.dirs]
        return str(variables) + " --> " + str(self.sym_func)


    def __mul__(self, other):
        if type(other) is not Function:
            other_func = Function(other, dirs=self.dirs, allow_sym=True)
            return self.__mul__(other_func)
        assert self.dirs == other.dirs
        new_sym_func = self.sym_func * other.sym_func
        allow_sym = self.allow_sym or other.allow_sym
        return Function(new_sym_func, dirs=self.dirs, allow_sym=allow_sym)

    def __add__(self, other):
        if type(other) is not Function:
            other_func = Function(other, dirs=self.dirs, allow_sym=True)
            return self.__add__(other_func)
        assert self.dirs == other.dirs
        new_sym_func = self.sym_func + other.sym_func
        allow_sym = self.allow_sym or other.allow_sym
        return Function(new_sym_func, dirs=self.dirs, allow_sym=allow_sym)

    __str__ = __repr__
    __rmul__ = __mul__
    __radd__ = __add__

    def as_xyz(self):
        result = self.sym_func
        for i in range(len(self.x_sub)):
            result = result.subs(self.x_sub[i], self.xyz[i])
        return result

    def ccode(self):
        function = sym.ccode(self.sym_func)
        for i in range(self.dim):
            function = re.sub(r'\bx{}'.format(self.dirs[i]),
                              'v[{}]'.format(i), function)
        return function

    def split(self):
        variables = [self.x_sub[d] for d in self.dirs]
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
