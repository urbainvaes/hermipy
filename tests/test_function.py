import hermite.function as func
import sympy as sym
import unittest

# import math
# import numpy as np
# import numpy.polynomial.hermite_e as herm
# import numpy.linalg as la

x, y, z = sym.symbols('x y z')

v = [sym.symbols('v[{}]'.format(i)) for i in range(3)]


class TestConstructor(unittest.TestCase):

    def test_construct_from_string(self):
        function_str1 = 'x*cos(y) + exp(z)**2'
        function_str2 = 'v0*cos(v1) + exp(v2)**2'
        function_sym1 = x*sym.cos(y) + sym.exp(z)**2
        function_sym2 = v[0]*sym.cos(v[1]) + sym.exp(v[2])**2
        function1 = func.Function(function_str1)
        function2 = func.Function(function_str2)
        function3 = func.Function(function_sym1)
        function4 = func.Function(function_sym2)
        self.assertTrue(function1 == function2)
        self.assertTrue(function2 == function3)
        self.assertTrue(function3 == function4)

    def test_as_string(self):
        function_str1 = 'v0*cos(v1) + exp(2*v2)'
        function_str2 = 'v[0]*cos(v[1]) + exp(2*v[2])'
        function_str3 = 'x*cos(y) + exp(2*z)'
        function = func.Function(function_str1)
        self.assertTrue(function.as_string(format='sub') == function_str1)
        self.assertTrue(function.as_string(format='array') == function_str2)
        self.assertTrue(function.as_string(format='xyz') == function_str3)
