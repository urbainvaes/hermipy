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

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "hermite/types.hpp"
#include "hermite/sparse.hpp"

namespace hermite {
namespace matrix
{
    inline void set(spmat & input, u_int i, u_int j, double val)
    {
        input(i, j) = val;
    }

    inline void set(mat & input, u_int i, u_int j, double val)
    {
        input[i][j] = val;
    }

    inline u_int size1(const mat & input)
    {
        return input.size();
    }

    inline u_int size2(const mat & input)
    {
        return input[0].size();
    }

    inline u_int size1(const spmat & input)
    {
        return input.size1();
    }

    inline u_int size2(const spmat & input)
    {
        return input.size2();
    }

    inline double get(const spmat & input, u_int i, u_int j)
    {
        return input(i, j);
    }

    inline double get(const mat & input, u_int i, u_int j)
    {
        return input[i][j];
    }


    template<typename T> T convert(const mat & input);
    template<typename T> T convert(const spmat & input);
    template<typename T> T construct(u_int size1, u_int size2);
    template<typename T> T eye(u_int size)
    {
        T result = construct<T> (size, size);
        for (u_int i = 0; i < size; i++)
        {
            set(result, i, i, 1.);
        }
        return result;
    }
}}

#endif
