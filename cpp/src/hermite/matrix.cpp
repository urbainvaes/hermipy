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

#include "hermite/matrix.hpp"

namespace hermite { namespace matrix {

    // Template specialization
    template <> mat construct(u_int size1, u_int size2)
    {
        return mat(size1, vec(size2, 0.));
    }

    template <> spmat construct(u_int size1, u_int size2)
    {
        return spmat(size1, size2);
    }

    template <> cmat construct(u_int size1, u_int size2)
    {
        auto dims = boost::extents[size1][size2];
        return cmat(dims);
    }

    template<typename T> T convert(const mat & input);
    template<typename T> T convert(const spmat & input);

    template <> mat convert (const mat & input)
    {
        return input;
    }

    template <> spmat convert (const spmat & input)
    {
        return input;
    }

    template <> mat convert(const spmat & input)
    {
        return full(input);
    }

    template <> spmat convert(const mat & input)
    {
        return to_spmat(input);
    }
}}
