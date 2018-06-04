import hermite.function as func
import sympy as sym
import unittest

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

    def test_auto_dim(self):
        newf = func.Function
        f1, f2, f3 = newf('x'), newf('x*y'), newf('x*y*z')
        self.assertTrue(f1.dim == 1 and f2.dim == 2 and f3.dim == 3)


class TestSplit(unittest.TestCase):

    def test_simple_split(self):
        newf = func.Function
        f = func.Function('x*y**2 + 5*exp(z)*y + 3*log(x*z)')
        split = f.split()
        self.assertTrue(f.dim == 3)
        self.assertTrue(len(split) == 3)
        for i in range(3):
            if split[i][-1] == 1:
                self.assertTrue(len(split[i]) == 4)
                self.assertTrue(split[i][0] == newf('x'))
                self.assertTrue(split[i][1] == newf('x*x'))
            if split[i][-1] == 5:
                self.assertTrue(len(split[i]) == 4)
                self.assertTrue(split[i][0] == newf('1'))
                self.assertTrue(split[i][1] == newf('x'))
                self.assertTrue(split[i][2] == newf('exp(x)'))
            if split[i][-1] == 3:
                self.assertTrue(len(split[i]) == 2)
                self.assertTrue(split[i][0] == newf('log(x*z)'))
