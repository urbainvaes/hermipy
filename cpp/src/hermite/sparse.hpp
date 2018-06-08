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

#ifndef SPARSE_H
#define SPARSE_H

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "hermite/types.hpp"

namespace bnu = boost::numeric::ublas;
namespace hermite {

mat row_col_val(const spmat & input);
mat full(const spmat & input);
spmat to_spmat(const mat & input);

spmat to_spmat(
        const vec & data,
        const ivec & indices,
        const ivec & indptr,
        u_int size1,
        u_int size2);

}

#endif
