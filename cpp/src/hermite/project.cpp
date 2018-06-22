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

#include <iostream>
#include "hermite/io.hpp"
#include "hermite/iterators.hpp"
#include "hermite/lib.hpp"
#include "hermite/matrix.hpp"
#include "hermite/project.hpp"
#include "hermite/templates.hpp"
#include "hermite/types.hpp"

using namespace std;

namespace hermite
{

template<typename Iterator>
vec project(const vec & input, u_int dim, const ivec & dirs)
{
    u_int n_polys = input.size();
    u_int dim_sub = dirs.size();
    u_int degree = Iterator::s_bissect_degree(dim, n_polys);
    u_int polys_sub = Iterator::s_size(dim_sub, degree);
    vec results(polys_sub, 0.);

    Iterator m(dim, degree),
      m_sub(dim_sub, degree);

    u_int i;
    for (i = 0, m.reset(); !m.isFull(); i++, m.increment())
    {
        if (isAligned(m.get(), dirs))
        {
            ivec sub = extract(m.get(), dirs);
            u_int ind = m_sub.index(sub);
            results[ind] = input[i];
        }
    }
    return results;
}

vec project_vec_nd(const vec & input, u_int dim, const ivec & dirs, std::string index_set)
{
    if (index_set == "cross")
    {
        return project<Cross_iterator>(input, dim, dirs);
    }
    else if (index_set == "triangle")
    {
        return project<Triangle_iterator>(input, dim, dirs);
    }
    else if (index_set == "cube")
    {
        return project<Cube_iterator>(input, dim, dirs);
    }
    else
    {
        std::cerr << "Invalid index set!" << std::endl;
        exit(1);
    }
}

vec project_vec_1d(const vec & input, u_int dim, u_int dir, std::string index_set)
{
    ivec dirs(1, dir);
    return project_vec_nd(input, dim, dirs, index_set);
}

template<typename Iterator, typename M>
M project(const M & input, u_int dim, const ivec & dirs)
{
    u_int n_polys = matrix::size1(input);
    u_int dim_sub = dirs.size();
    u_int degree = Iterator::s_bissect_degree(dim, n_polys);
    u_int polys_sub = Iterator::s_size(dim_sub, degree);
    M results = matrix::construct<M>(polys_sub, polys_sub);

    Iterator m1(dim, degree), m2(dim, degree),
      m_sub(dim_sub, degree);

    u_int i,j;
    for (i = 0, m1.reset(); !m1.isFull(); i++, m1.increment())
    {
        if (isAligned(m1.get(), dirs))
        {
            ivec sub1 = extract(m1.get(), dirs);
            u_int ind1 = m_sub.index(sub1);
            for (j = 0, m2.reset(); !m2.isFull(); j++, m2.increment())
            {
                if (isAligned(m2.get(), dirs))
                {
                    ivec sub2 = extract(m2.get(), dirs);
                    u_int ind2 = m_sub.index(sub2);
                    matrix::set(results, ind1, ind2, matrix::get(input, i, j));
                }
            }
        }
    }
    return results;
}

template<typename M>
M project_mat_nd(const M & input, u_int dim, const ivec & dirs, std::string index_set)
{
    if (index_set == "cross")
    {
        return project<Cross_iterator,M>(input, dim, dirs);
    }
    else if (index_set == "triangle")
    {
        return project<Triangle_iterator,M>(input, dim, dirs);
    }
    else if (index_set == "cube")
    {
        return project<Cube_iterator,M>(input, dim, dirs);
    }
    else
    {
        std::cerr << "Invalid index set!" << std::endl;
        exit(1);
    }
}

template <typename M>
M project_mat_1d(const M & input, u_int dim, u_int dir, std::string index_set)
{
    ivec dirs(1, dir);
    return project_mat_nd(input, dim, dirs, index_set);
}

template mat project_mat_1d(const mat & input, u_int dim, u_int dir, std::string index_set);
template spmat project_mat_1d(const spmat & input, u_int dim, u_int dir, std::string index_set);
template mat project_mat_nd(const mat & input, u_int dim, const ivec & dirs, std::string index_set);
template spmat project_mat_nd(const spmat & input, u_int dim, const ivec & dirs, std::string index_set);

}
