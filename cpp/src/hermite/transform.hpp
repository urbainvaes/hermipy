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

#ifndef HERMITE_TRANSFORM_H
#define HERMITE_TRANSFORM_H

#include <string>
#include "hermite/types.hpp"

namespace hermite {

vec transform(
        u_int degree,
        vec const & f_grid,
        mat const & nodes,
        mat const & weights,
        bool forward,
        std::string index_set = "triangle");

double integrate(
        vec const & f_grid,
        mat const & nodes,
        mat const & weights);
}

#endif
