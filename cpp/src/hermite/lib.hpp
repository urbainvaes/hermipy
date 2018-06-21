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

#ifndef LIB_H
#define LIB_H

#include <string>
#include <functional>
#include "hermite/types.hpp"

namespace hermite
{

// Check if multi-index is aligned to another
bool isAligned(const ivec & m, u_int dir);
bool isAligned(const ivec & m, const ivec & dirs);


// Extract sub-vector
template<typename T>
std::vector<T> extract (const std::vector<T> & input, const ivec & indices)
{
    std::vector<T> result(indices.size(), 0);
    for (u_int i = 0; i < indices.size(); i++)
    {
        result[i] = input[indices[i]];
    }
    return result;
}

u_int hash_multi_ind(const ivec & v, int degree);
std::string hash_print(const ivec & v);
u_int pos_bissect(u_int image, std::function<u_int(u_int)>, u_int max);

}

#endif
