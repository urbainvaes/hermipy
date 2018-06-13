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

// FIXME: Last multi-index twice? (urbain, 13 Jun 2018)


#include <unordered_map>
#include "hermite/io.hpp"
#include "hermite/types.hpp"
#include "hermite/iterators.hpp"

#include <boost/math/special_functions/binomial.hpp>

namespace hermite
{

u_int hash_multi_ind(ivec v, int degree)
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

#define MAX(i, j) (i > j ? i : j)
bool Hyperbolic_cross_increment(ivec & multi_index, u_int dim, u_int upper_bound)
{
    u_int i;

    for (i = 0; i < dim; i++)
    {
        if (multi_index[i] == 1)
            multi_index[i] = 0;

        else if (multi_index[i] == 0)
        {
            multi_index[i] = 1;
            return false;
        }
    }

    for (i = 0; i < dim - 1; i++)
    {
        if (multi_index[i + 1] == 0)
            continue;

        bool found_factor = false;
        u_int factor;
        for (u_int f : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31})
        {
            if (multi_index[i + 1] % f == 0)
            {
                found_factor = true;
                factor = f;
                break;
            }
        }
        if (!found_factor)
            factor = multi_index[i+1];

        bool is_prime = factor == multi_index[i+1];

        multi_index[i] = MAX(1, multi_index[0]) * factor;
        multi_index[i + 1] /= factor;

        if (is_prime)
            multi_index[i + 1] = 0;

        break;
    }

    if (i == dim - 1)
    {
        if (multi_index[0] == upper_bound)
            return true;
        multi_index[i] = 1 + MAX(1, multi_index[0]);
    }

    if (i > 0)
        multi_index[0] = 0;

    return false;
}

Hyperbolic_cross_iterator::Hyperbolic_cross_iterator(u_int dim, u_int upper_bound)
    : Vector_iterator(dim), upper_bound(upper_bound), index_list(0)
{
    ivec m(dim, 0);
    u_int ind = 0;
    do
    {
        list.push_back(m);
        u_int hash = hash_multi_ind(m, upper_bound);
        hash_table.insert(std::pair<u_int, u_int>(hash, ind++));
    }
    while (!Hyperbolic_cross_increment(m, dim, upper_bound));
}

void Hyperbolic_cross_iterator::increment()
{
    if (index_list == list.size())
    {
        full = true;
        return;
    }

    multi_index = list[++index_list];
}

void Multi_index_iterator::increment()
{
    u_int i;
    for (i = 0; i < dim - 1; i++)
    {
        if (multi_index[i + 1]  != 0)
        {
            multi_index[i + 1] -= 1;
            break;
        }
    }

    if (i == dim - 1 && multi_index[0] == upper_bound)
    {
        full = true;
        return;
    }

    multi_index[i] = 1 + multi_index[0];

    if (i > 0)
    {
        multi_index[0] = 0;
    }
}

void Hyper_cube_iterator::increment()
{
    unsigned int i = dim - 1;
    while(multi_index[i] == upper_bounds[i] - 1 && i > 0)
    {
        multi_index[i] = 0;
        i -= 1;
    }
    if(i == 0 && multi_index[0] == upper_bounds[0] - 1)
    {
        full = true;
    }
    else
    {
        multi_index[i] += 1;
    }
}

imat Hyper_cube_iterator::list(const ivec & upper_bounds)
{
    imat result;
    for(Hyper_cube_iterator m(upper_bounds); !m.isFull(); m.increment())
    {
        result.push_back(m.get());
    }
    return result;
}

imat Multi_index_iterator::list(u_int dim, u_int upper_bound)
{
    imat result;
    for(Multi_index_iterator m(dim, upper_bound); !m.isFull(); m.increment())
    {
        result.push_back(m.get());
    }
    return result;
}

u_int Multi_index_iterator::size(u_int degree, u_int dim)
{
    using boost::math::binomial_coefficient;
    return (u_int) binomial_coefficient<double> (degree + dim, dim);
}

u_int Multi_index_iterator::index(const ivec & m_vec)
{
    using boost::math::binomial_coefficient;
    u_int sum = 0, result = 0;
    for (u_int i = 0; i < m_vec.size(); i++)
    {
        sum += m_vec[i];
        if (sum > 0)
        {
            result += (u_int) binomial_coefficient<double> (i + sum, i + 1);
        }
    }
    return result;
}

u_int Hyperbolic_cross_iterator::index(const ivec & m_vec)
{
    u_int hash_mvec = hash_multi_ind(m_vec, upper_bound);
    return hash_table.at(hash_mvec);
}

u_int Hyper_cube_iterator::index(const ivec & m_vec)
{
    u_int result = 0,
          factor = 1;

    for (u_int i = 0; i < dim; i++)
    {
        u_int position = dim - 1 - i;
        result += factor * m_vec[position];
        factor *= upper_bounds[position];
    }
    return result;
}

}
