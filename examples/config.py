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

from hermipy.equations import McKean_Vlasov as equation
import sympy as sym

# Configuration dicionaries
misc, eq, num = {}, {}, {}

# Variables and function
x, y, f = equation.x, equation.y, equation.f

# Short-hand notation
r = sym.Rational

# Configuration of numerical method
num['degree'] = 100  # degree of approximation
num['n_points_num'] = 2*num['degree'] + 1  # (*2 for varf)
num['μx'] = r(1, 5)
num['σx'] = r(1, 10)
num['λ'] = r(1, 2)

# Scalar parameters of the equation
eq['β'] = r(2)
eq['ε'] = r(1)
eq['γ'] = r(0)
eq['θ'] = r(0)

# Functional parameters of the equation
# eq['Vp'] = x**4/4 - x**2/2

sx = sym.symbols('sx', real=True, positive=True)
eq['sx'] = r(1)
eq['Vp'] = r(.5)*x*x/sx

# Miscellaneous parameters
misc['cache'] = False
misc['debug'] = False
misc['parallel'] = False
misc['tensorize'] = True
misc['trails'] = True
misc['plots'] = True
misc['symbolic'] = 0  # Values 0, 1, 2
