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

#ifndef PROJECT_H
#define PROJECT_H

#include "hermite/types.hpp"

namespace hermite
{

vec project_vec_1d(const vec & input, u_int dim, u_int dir, std::string index_set);
vec project_vec_nd(const vec & input, u_int dim, const ivec & dirs, std::string index_set);

template <typename M> M project_mat_1d(const M & input, u_int dim, u_int dir, std::string index_set);
template <typename M> M project_mat_nd(const M & input, u_int dim, const ivec & dirs, std::string index_set);

}

#endif
