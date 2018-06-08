/*
 * Copyright (C) 2018 Urbain Vaes
 *
 * This file is part of hermipy, a python/C++ library for automating the
 * Hermite Galerkin method.
 *
 * hermipy is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hermipy is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>

#include "hermite/hermite.hpp"

using namespace std;

namespace hermite {

void hermite_eval(double x,
        u_int degree,
        vec & values)
{

    values[0] = 1;
    values[1] = x;

    for (unsigned i = 1; i < degree; i++)
    {
        values[i+1] = REC_A(i)*x*values[i] - REC_B(i)*values[i-1];
    }
}}
