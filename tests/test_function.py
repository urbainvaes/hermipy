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
        function_str = 'x*cos(x) + exp(x)**2'
        function_sym = v[0]*sym.cos(v[0]) + sym.exp(v[0])**2
        function = func.Function(function_str)
        print(function.sym_func)
        self.assertTrue(function.sym_func == function_sym)
