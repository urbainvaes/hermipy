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

import sympy
import re


class Function():

    # xyz notation
    xyz = list(sympy.symbols('x y z v w', real=True))

    # subscript notation
    x_sub = [sympy.symbols('x{}'.format(i), real=True)
             for i in range(len(xyz))]

    conv = {}
    for i in range(len(xyz)):
        conv[str(xyz[i])] = x_sub[i]
        conv[str(x_sub[i])] = x_sub[i]

    def __init__(self, expr, dirs=None, dim=None):

        if isinstance(expr, Function):
            self.sym = expr.sym
            expr.dirs = expr.dirs
            return

        if not isinstance(expr, tuple(sympy.core.all_classes)):
            expr = sympy.sympify(expr, locals=self.conv)

        for s in expr.free_symbols:
            if str(s) in self.conv:
                expr = expr.subs(s, self.conv[str(s)])

        self.sym = expr
        variables = expr.free_symbols.intersection(self.x_sub)
        if dirs is not None:
            self.dirs = dirs
            assert dim is None
            assert self.dirs == sorted(self.dirs)
        elif dim is not None:
            self.dirs = list(range(dim))
        else:
            # We do not minimize the number of dimensions
            dim = 0 if variables == set() else \
                  max(self.x_sub.index(s) for s in variables) + 1
            self.dirs = list(range(dim))
        assert variables.issubset({self.x_sub[d] for d in self.dirs})

    def __eq__(self, other):
        return self.dirs == other.dirs and \
            self.sym == other.sym

    def __repr__(self):
        variables = [self.x_sub[d] for d in self.dirs]
        return str(variables) + " --> " + str(self.sym)

    def __mul__(self, other):
        if type(other) is not Function:
            other = Function(other, dirs=self.dirs)
        assert self.dirs == other.dirs
        sym = self.sym * other.sym
        return Function(sym, dirs=self.dirs)

    def __add__(self, other):
        if type(other) is not Function:
            other = Function(other, dirs=self.dirs)
        assert self.dirs == other.dirs
        sym = self.sym + other.sym
        return Function(sym, dirs=self.dirs)

    __str__ = __repr__
    __rmul__ = __mul__
    __radd__ = __add__

    def as_xyz(self):
        result = self.sym
        for i in range(len(self.x_sub)):
            result = result.subs(self.x_sub[i], self.xyz[i])
        return result

    def ccode(self):
        function = sympy.ccode(self.sym)
        for i in range(len(self.dirs)):
            function = re.sub(r'\bx{}'.format(self.dirs[i]),
                              'v[{}]'.format(i), function)
        return function

    def split(self):
        variables = [self.x_sub[d] for d in self.dirs]
        is_add = isinstance(self.sym, sympy.add.Add)
        add_terms = self.sym.args if is_add else [self.sym]
        to_return = []
        for term in add_terms:
            is_mul = isinstance(term, sympy.mul.Mul)
            mul_terms = term.args if is_mul else [term]
            result = [sympy.Rational('1')] * (len(self.dirs) + 1)
            tensorizable = True
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
                for i in range(len(self.dirs)):
                    result[i] = Function(result[i], dirs=[self.dirs[i]])
            else:
                func = Function(term/result[-1], dirs=self.dirs)
                result = [func, result[-1]]
            to_return.append(result)
        return to_return

    @staticmethod
    def sanitize(expr, max_denom=1e8):
        result, expr = sympy.Integer(0), (sympy.Integer(0) + expr).expand()
        is_add = isinstance(expr, sympy.add.Add)
        add_terms = expr.args if is_add else [expr]
        for aterm in add_terms:
            result_tmp = sympy.Integer(1)
            is_mul = isinstance(aterm, sympy.mul.Mul)
            mul_terms = aterm.args if is_mul else [aterm]
            for mterm in mul_terms:
                if type(mterm) is sympy.numbers.Float:
                    if abs(mterm) < 1e-14:
                        mterm = 0
                    mterm = sympy.Rational(mterm).limit_denominator(max_denom)
                result_tmp *= mterm
            result += result_tmp
        return result
