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

import hermipy.function as func
import sympy
import itertools
import math


class Operator():

    f = sympy.Function('f')

    def __init__(self, expr):

        if isinstance(expr, Operator):
            self.sym = expr.sym
            self.dirs = expr.dirs
            return

        function = func.Function(expr)
        self.sym = function.sym
        variables = self.sym.free_symbols.intersection(func.Function.x_sub)
        self.dirs = sorted([func.Function.x_sub.index(s) for s in variables])

        if not self.sym.has(self.f):
            raise TypeError("Argument is not an operator on f" + str(self.sym))

    def __eq__(self, other):
        return self.dirs == other.dirs and \
            self.sym == other.sym

    def __repr__(self):
        variables = [func.Function.x_sub[d] for d in self.dirs]
        return str(self.f(*variables)) + " --> " + str(self.sym)

    def __mul__(self, other):
        if type(other) is not func.Function:
            other = func.Function(other, dirs=self.dirs)
        assert self.dirs == other.dirs
        variables = [func.Function.x_sub[d] for d in self.dirs]
        sym = self.sym.subs(self.f(*variables), other.sym)
        return Operator(sym)

    def __call__(self, arg):

        if not isinstance(arg, (func.Function, Operator)):
            try:
                arg = Operator(arg)
            except TypeError:
                arg = func.Function(arg)

        variables = [func.Function.x_sub[d] for d in self.dirs]
        sym = self.sym.subs(self.f(*variables), arg.sym)
        return Operator(sym) if type(arg) is Operator \
            else func.Function(sym, dirs=self.dirs)

    def __add__(self, other):
        if type(other) is not Operator:
            other = Operator(other, dirs=self.dirs)
        assert self.dirs == other.dirs
        sym = self.sym + other.sym
        return Operator(sym, dirs=self.dirs)

    __str__ = __repr__
    __radd__ = __add__
    __rmul__ = __mul__

    def as_xyz(self):
        return self.sym.subs(((func.Function.x_sub[i], func.Function.xyz[i])
                              for i in range(len(func.Function.x_sub))))

    def split(self):
        MAX_ORDER = 5
        variables = [func.Function.x_sub[d] for d in self.dirs]
        unknown, dim = self.f(*variables), len(variables)
        result, rem = [], self.sym.expand()
        mult = list(m for m in itertools.product(range(MAX_ORDER + 1),
                    repeat=dim) if sum(m) <= MAX_ORDER)
        for m in mult:
            if rem == 0:
                result.append(0)
                continue
            test, der = 1, unknown
            for i, v in zip(m, variables):
                test *= v**i/math.factorial(i)
                der = sympy.diff(der, v, i)
            remargs = rem.args if isinstance(rem, sympy.add.Add) else [rem]
            term, rem = 0, 0
            for arg in remargs:  # Convoluted to avoid rounding errors
                termarg = arg.subs(unknown, test).doit()
                if termarg == 0:
                    rem += arg
                else:
                    term += termarg
            if isinstance(term, tuple(sympy.core.all_classes)):
                term = sympy.simplify(term)
            result.append(term)
        assert rem == 0
        result = []
        nz =  [(m, r) for m, r in zip(mult, result) if r != 0]
        for m, r in nz:
            split_coeff = func.Function(nz[i], dirs=self.dirs).split()
            tensorizable = len(split_coeff) == len(self.dirs) + 1
            if tensorizable:
                split_op = []
                for i in range(dim):
                    derivative = unknown.diff(*([variables[i]]*m[i]))
                    split_op.append(split_coeff[i]*derivative)
                split_op.append(split_coeff[-1])
