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
import sympy as sym
import unittest

x, y, z = sym.symbols('x y z')

v = [sym.symbols('x{}'.format(i)) for i in range(3)]


class TestFunction(unittest.TestCase):

    def test_construct_from_string(self):
        function_str1 = 'x*cos(y) + exp(z)**2'
        function_str2 = 'x0*cos(x1) + exp(x2)**2'
        function_sym1 = x*sym.cos(y) + sym.exp(z)**2
        function_sym2 = v[0]*sym.cos(v[1]) + sym.exp(v[2])**2
        function1 = func.Function(function_str1)
        function2 = func.Function(function_str2)
        function3 = func.Function(function_sym1)
        function4 = func.Function(function_sym2)
        self.assertTrue(function1 == function2)
        self.assertTrue(function2 == function3)
        self.assertTrue(function3 == function4)

    def test_as_xyz(self):
        function_str1 = 'x0*cos(x1) + exp(2*x2)'
        function_str2 = 'x*cos(y) + exp(2*z)'
        function = func.Function(function_str1)
        self.assertTrue(str(function.sym) == function_str1)
        self.assertTrue(str(function.as_xyz()) == function_str2)

    def test_auto_dim(self):
        newf = func.Function
        f1, f2, f3 = newf('x'), newf('x*y'), newf('x*y*z')
        self.assertTrue(len(f1.dirs) == 1 and len(f2.dirs) == 2
                        and len(f3.dirs) == 3)

    def test_simple_split(self):
        newf = func.Function
        f = func.Function('x*y**2 + 5*exp(z)*y + 3*log(x*z)')
        split = f.split()
        self.assertTrue(len(f.dirs) == 3)
        self.assertTrue(len(split) == 3)
        for i in range(3):
            if split[i][-1] == 1:
                self.assertTrue(len(split[i]) == 4)
                self.assertTrue(split[i][0] == newf('x', dirs=[0]))
                self.assertTrue(split[i][1] == newf('y*y', dirs=[1]))
            if split[i][-1] == 5:
                self.assertTrue(len(split[i]) == 4)
                self.assertTrue(split[i][0] == newf('1', dirs=[0]))
                self.assertTrue(split[i][1] == newf('y', dirs=[1]))
                self.assertTrue(split[i][2] == newf('exp(z)', dirs=[2]))
            if split[i][-1] == 3:
                self.assertTrue(len(split[i]) == 2)
                self.assertTrue(split[i][0] == newf('log(x*z)', dim=3))


if __name__ == '__main__':
    unittest.main()
