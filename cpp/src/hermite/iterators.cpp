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

#include <assert.h>
#include <unordered_map>
#include "hermite/io.hpp"
#include "hermite/types.hpp"
#include "hermite/iterators.hpp"
#include <boost/math/special_functions/binomial.hpp>

#define MAX(i, j) (i > j ? i : j)
#define DEGREE_BISSECT_MAX 150

namespace hermite
{
// Grid iterators (used for N-dim for loops)
// Hypercube iterator {{{
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

// }}}

// Multi-index iterators
// Parent class {{{
u_int Multi_index_iterator::index(const ivec & m_vec)
{
    std::string hash_mvec = hash_print(m_vec);
    if (hash_table.find(hash_mvec) == hash_table.end())
    {
        std::cout << "Element not found in hash table!" << std::endl;
        std::cout << "--> " << m_vec << std::endl;
        std::cout << "--> " << list << std::endl;
    }
    return hash_table.at(hash_mvec);
}

void Multi_index_iterator::increment()
{
    if (index_list == list.size() - 1)
    {
        full = true;
        return;
    }

    multi_index = list[++index_list];
}
// }}}
// Triangle iterator {{{
bool Triangle_iterator::s_increment(ivec & multi_index, u_int dim, u_int upper_bound)
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
        return true;

    multi_index[i] = 1 + multi_index[0];

    if (i > 0)
    {
        multi_index[0] = 0;
    }

    return false;
}

Triangle_iterator::Triangle_iterator(u_int dim, u_int upper_bound)
    : Multi_index_iterator(dim), upper_bound(upper_bound)
{
    #ifdef DEBUG
    std::cout << "Creating triangle iterator with dim = "
              << dim << ", upper_bound = " << upper_bound << std::endl;
    #endif
    ivec m(dim, 0);
    u_int ind = 0;
    do
    {
        list.push_back(m);
        std::string hash = hash_print(m);
        hash_table.insert(std::pair<std::string, u_int>(hash, ind++));
    }
    while (!Triangle_iterator::s_increment(m, dim, upper_bound));
}

u_int Triangle_iterator::s_index(const ivec & m_vec)
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

u_int Triangle_iterator::s_size(u_int degree, u_int dim)
{
    using boost::math::binomial_coefficient;
    return (u_int) binomial_coefficient<double> (degree + dim, dim);
}

imat Triangle_iterator::s_list(u_int dim, u_int upper_bound)
{
    imat result;
    for(Triangle_iterator m(dim, upper_bound); !m.isFull(); m.increment())
    {
        result.push_back(m.get());
    }
    return result;
}

u_int Triangle_iterator::s_bissect_degree(u_int dim, u_int n_polys)
{
    u_int degree_1 = 0, degree_2 = DEGREE_BISSECT_MAX;

    int img_1 = (int) s_size(degree_1, dim) - (int) n_polys;
    int img_2 = (int) s_size(degree_2, dim) - (int) n_polys;

    if (img_1 > 0 || img_2 < 0)
    {
        std::cout << "Can't find degree, Invalid arguments!" << std::endl;
        exit(0);
    }

    if (img_1 == 0)
        return degree_1;

    if (img_2 == 0)
        return degree_2;

    while (true)
    {
        u_int new_degree = (degree_1 + degree_2)/2;
        int new_img = (int) s_size(new_degree, dim) - (int) n_polys;

        if (new_img < 0)
        {
            degree_1 = new_degree;
            img_1 = new_img;
        }
        else if (new_img > 0)
        {
            degree_2 = new_degree;
            img_2 = new_img;
        }
        else
        {
            return new_degree;
        }
    }
}

u_int Triangle_iterator::s_find_dim(u_int degree, u_int n_polys)
{
    using boost::math::binomial_coefficient;

    u_int dim = 0;
    u_int n_dim = (u_int) binomial_coefficient<double> (degree + dim, dim);

    while (n_dim < n_polys)
    {
        dim += 1;
        n_dim = (u_int) binomial_coefficient<double> (degree + dim, dim);
    }

    if (n_dim != n_polys)
    {
        std::cerr << "Dimension not found, exiting..." << std::endl;
        exit(1);
    }

    return dim;
}

// }}}
// Hyperbolic cross iterator {{{
bool Cross_iterator::s_increment(ivec & multi_index, u_int dim, u_int upper_bound)
{
    assert(upper_bound > 0);
    u_int i;

    // FIXME: Bug when degree == 1
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

    if (upper_bound == 1)
    {
        return true;
    }

    for (i = 0; i < dim - 1; i++)
    {
        if (multi_index[i + 1] == 0)
            continue;

        u_int factor = 0;
        u_int product = MAX(1, multi_index[0]) * multi_index[i + 1];
        for (u_int f = MAX(1, multi_index[0]) + 1; f <= product; f++)
        {
            if (product % f == 0)
            {
                factor = f;
                break;
            }
        }
        multi_index[i] = factor;
        multi_index[i + 1] = factor == product ? 0 : product / factor;
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

Cross_iterator::Cross_iterator(u_int dim, u_int upper_bound)
    : Multi_index_iterator(dim), upper_bound(upper_bound)
{
    #ifdef DEBUG
    std::cout << "Creating cross iterator with dim = "
              << dim << ", upper_bound = " << upper_bound << std::endl;
    #endif

    ivec m(dim, 0);
    u_int ind = 0;
    do
    {
        list.push_back(m);
        std::string hash = hash_print(m);
        hash_table.insert(std::pair<std::string, u_int>(hash, ind++));
    }
    while (!Cross_iterator::s_increment(m, dim, upper_bound));
}

imat Cross_iterator::s_list(u_int dim, u_int upper_bound)
{
    imat result;
    for(Cross_iterator m(dim, upper_bound); !m.isFull(); m.increment())
    {
        result.push_back(m.get());
    }
    return result;
}

u_int Cross_iterator::s_bissect_degree(u_int dim, u_int n_polys)
{
    Cross_iterator m(dim, DEGREE_BISSECT_MAX);
    assert (n_polys <= m.list.size());
    return m.list[n_polys - 1][0];
}

u_int Cross_iterator::s_size(u_int degree, u_int dim)
{
    return Cross_iterator(dim, degree).size();
}
// }}}
}
