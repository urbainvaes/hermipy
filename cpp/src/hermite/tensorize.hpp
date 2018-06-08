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

#ifndef TENSORIZE_H
#define TENSORIZE_H

#include "hermite/types.hpp"
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace hermite {

// Tensorize vector(s)
vec tensorize(const vec & input, u_int dim, u_int dir);
vec tensorize(const mat &);
vec tensorize(const mat & inputs, const imat & dirs);

// Tensorize matrix
template<typename T, typename M>
T tensorize(const M & input, u_int dim, u_int dir);

template <typename T, typename M>
T tensorize(const std::vector<M> & inputs);

template <typename T, typename S>
T tensorize(const std::vector<S> & inputs, const imat & dirs);

}

#endif
