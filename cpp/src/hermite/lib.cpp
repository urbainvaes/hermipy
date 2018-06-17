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

#include "hermite/iterators.hpp"
#include "hermite/lib.hpp"
#include "hermite/types.hpp"

using namespace std;

namespace hermite
{

bool isAligned(const ivec & m, u_int dir)
{
    for (u_int j = 0; j < m.size(); j++)
    {
        if (m[j] != 0 && j != dir)
        {
            return false;
        }
    }
    return true;
}

bool isAligned(const ivec & m, const ivec & dirs)
{
    for (u_int j = 0; j < m.size(); j++)
    {
        bool in_dirs = false;
        for(u_int k = 0; k < dirs.size(); k++)
        {
            if (j == dirs[k])
            {
                in_dirs = true;
            }
        }

        if (m[j] != 0 && !in_dirs)
        {
            return false;
        }
    }
    return true;
}

u_int hash_multi_ind(const ivec & v, int degree)
{
    u_int base = degree + 1;
    u_int result = 0;
    u_int unit = 1;
    for(u_int i = 0; i < v.size(); i++)
    {
        result += v[i]*unit;
        unit *= base;
    }
    return result;
}

string hash_print(const ivec & v)
{
    string hash = "";
    for (u_int i = 0; i < v.size(); ++i)
    {
        hash += std::to_string(v[i]) + "-";
    }
    return hash;
}

}
