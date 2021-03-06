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

# Classes in global namespace

from hermipy.quad import *
from hermipy.series import *
from hermipy.varf import *
from hermipy.function import *
from hermipy.operator import *
from hermipy.position import *

# Functions in global namespace
from hermipy.cache import cache

# Settings of the library
from hermipy.settings import *

# For convenience
x, y, z = Function.xyz[0:3]
f = Operator.f
