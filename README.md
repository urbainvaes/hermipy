# Hermipy

[![Documentation Status](https://readthedocs.org/projects/hermipy/badge/?version=latest)](https://hermipy.readthedocs.io/en/latest/?badge=latest)
[![Build Status](http://vaes.uk:8090/buildStatus/icon?job=hermipy)](http://vaes.uk:8090/job/hermipy/lastBuild/consoleText)
[![](http://vaes.uk:8090/job/hermipy/ws/coverage.svg)](http://vaes.uk:8090/job/hermipy/ws/coverage.txt)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a1588a65d35340959ecb5ca500e14121)](https://app.codacy.com/app/urbainvaes/hermipy?utm_source=github.com&utm_medium=referral&utm_content=urbainvaes/hermipy&utm_campaign=Badge_Grade_Dashboard)

Hermipy is a thin Python library that allows automating most of the operations involved in implementing a Hermite spectral method.
The library uses the data structures provided by the [NumPy](https://en.wikipedia.org/wiki/NumPy) library for linear algebra,
and it also depends on [SymPy](https://en.wikipedia.org/wiki/SymPy) for symbolic calculations.
Computationally intensive parts of the computations,
such as Hermite transforms,
are handed over to a C++ compiled component,
and [Boost](https://en.wikipedia.org/wiki/Boost_(C%2B%2B_libraries)) is used for efficient exchange of data with Python.
For more information, read the [documentation](https://hermipy.readthedocs.io/en/latest/).

# Todo

- Use third-party library to compute FFTs;

- Implement composite quadrature;

- Clean up data conversion between C++ and Python;

- Inherit sympy.function / sympy.operator instead of custom classes.

# License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
