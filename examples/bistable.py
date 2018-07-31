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
import hermipy

# Configuration dicionaries
misc, eq, num = {}, {}, {}

# Variables and function
x, y, f = equation.x, equation.y, equation.f

# Short-hand notation
r = sym.Rational

# Configuration of numerical method
num['degree'] = 40  # degree of approximation
num['n_points_num'] = 2*num['degree'] + 1  # (*2 for varf)
num['μx'] = r(0)
num['μy'] = r(0, 5)
num['σx'] = r(1, 25)
num['σy'] = r(1, 25)
num['λ'] = r(1, 2)
num['index_set'] = 'triangle'

# Scalar parameters of the equation
eq['β'] = r(8)
eq['ε'] = r(1, 2**6)
eq['γ'] = r(0)
eq['θ'] = r(0)

# Functional parameters of the equation
eq['Vp'] = x**4/4 - x**2/2
# eq['Vp'] = x**2/2
# eq['Vp'] = x**4

# Mean-zero
Z, m = 6.301119049538182, 0.8852269357209047
eq['Vy'] = (y-m)**4/4 - (y-m)**2/2 + (y-m)
# eq['Vy'] = y**4/4 - y**2/2
# eq['Vy'] = y**2/2

# rx, ry, d = sym.exp(-eq['Vp']), sym.exp(-eq['Vy']), num['degree']
# qy = hermipy.Quad.gauss_hermite(200, cov=[[.05]], dirs=[1])
# assert abs(qy.integrate(y*sym.exp(-eq['Vy']), flat=True)) < 1e-10

# qx.norm(qx.transform(rx, degree=d).eval(), qx.discretize(yx))
# qy.plot(sym.exp(-eq['Vy']))

# Miscellaneous parameters
misc['cache'] = True
misc['parallel'] = False
misc['tensorize'] = True
misc['trails'] = False
misc['verbose'] = False
misc['symbolic'] = 0  # Values 0, 1, 2
misc['plots'] = False
