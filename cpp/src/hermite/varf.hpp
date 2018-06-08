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

#ifndef HERMITE_VARF_H
#define HERMITE_VARF_H

#include "hermite/types.hpp"
#include "boost/multi_array.hpp"

namespace hermite {

cube triple_products_1d(int degree);

template<typename T>
T varf(
        u_int degree,
        vec const & input,
        mat const & nodes,
        mat const & weights);

mat varfd(
        u_int dim,
        u_int degree,
        u_int direction,
        const mat & var);

spmat varfd(
        u_int dim,
        u_int degree,
        u_int direction,
        const spmat & var);

}

#endif
