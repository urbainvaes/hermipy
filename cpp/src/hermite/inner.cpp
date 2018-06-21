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

#ifdef DEBUG
#include <iostream>
#endif

#include "hermite/inner.hpp"
#include "hermite/types.hpp"
#include "hermite/iterators.hpp"
#include "hermite/lib.hpp"

namespace hermite
{

// Precondition: Dirs in increasing order
template<typename Iterator>
vec inner(const vec & s1,
          const vec & s2,
          const ivec & dirs1,
          const ivec & dirs2)
{
    #ifdef DEBUG
    std::cout << "Entering in inner()." << std::endl;
    #endif

    u_int dim1 = dirs1.size(),
          dim2 = dirs2.size();

    u_int i = 0,
          j = 0;

    ivec dirs_inner,
         dirs_result;

    ivec inner_ind1,
         inner_ind2;

    std::vector<bool> in_dirs1;
    ivec free_ind;

    #ifdef DEBUG
    std::cout << "--> Starting while loop" << std::endl;
    #endif

    while (i < dim1 || j < dim2)
    {
        if (i == dim1 || (j < dim2 && dirs1[i] > dirs2[j]) )
        {
            in_dirs1.push_back(false);
            free_ind.push_back(j);
            dirs_result.push_back(dirs2[j]);
            j++;
        }
        else if (j == dim2 || (i < dim1 && dirs1[i] < dirs2[j]) )
        {
            in_dirs1.push_back(true);
            free_ind.push_back(i);
            dirs_result.push_back(dirs1[i]);
            i++;
        }
        else
        {
            assert (dirs1[i] == dirs2[j]);
            dirs_inner.push_back(dirs1[i]);
            inner_ind1.push_back(i); i++;
            inner_ind2.push_back(j); j++;
        }
    }

    u_int dim_result = dirs_result.size(),
          dim_inner = dirs_inner.size();

    #ifdef DEBUG
    std::cout << "Calculating degree with dim = " << dim1 << "," << dim2
        << " and n_polys = " << s1.size() << "," << s2.size() <<  std::endl;
    #endif

    u_int degree = Iterator::s_bissect_degree(dim1, s1.size());
    assert( degree == Iterator::s_bissect_degree(dim2, s2.size()) );

    #ifdef DEBUG
    std::cout << "Initializing iterators" <<  std::endl;
    #endif

    vec result(Iterator::s_size(degree, dim_result), 0.);

    Iterator m_result(dim_result, degree);
    Iterator m_inner(dim_inner, degree);
    Iterator it_m1(dim1, degree);
    Iterator it_m2(dim2, degree);

    #ifdef DEBUG
    std::cout << "--> Starting outer for loop" << std::endl;
    #endif

    for (i = 0; !m_result.isFull(); m_result.increment(), i++)
    {
        #ifdef DEBUG
        std::cout << "----> Value of m_result: " << m_result.get() << std::endl;
        #endif

        ivec m1(dim1),
             m2(dim2);

        for (u_int k = 0; k < dim_result; k++)
        {
            if (in_dirs1[k])
                m1[free_ind[k]] = m_result[k];
            else
                m2[free_ind[k]] = m_result[k];
        }

        #ifdef DEBUG
        std::cout << "----> Starting inner for loop." << std::endl;
        #endif
        for (j = 0, m_inner.reset(); !m_inner.isFull(); m_inner.increment(), j++)
        {
            #ifdef DEBUG
            std::cout << "------> Value of m_inner: " << m_inner.get() << std::endl;
            #endif

            for (u_int k = 0; k < dim_inner; k++)
            {
                m1[inner_ind1[k]] = m_inner[k];
                m2[inner_ind2[k]] = m_inner[k];
            }

            if (!it_m1.has(m1) || !it_m2.has(m2))
                continue;

            int ind1 = it_m1.index(m1),
                ind2 = it_m2.index(m2);

            result[i] += s1[ind1] * s2[ind2];
        }
    }

    #ifdef DEBUG
    std::cout << "Exiting inner()." << std::endl;
    #endif

    return result;
}

vec inner(const vec & s1,
          const vec & s2,
          const ivec & dirs1,
          const ivec & dirs2,
          std::string index_set)
{
    if (index_set == "cross")
    {
        return inner<Cross_iterator>(s1, s2, dirs1, dirs2);
    }
    else if (index_set == "triangle")
    {
        return inner<Triangle_iterator>(s1, s2, dirs1, dirs2);
    }
    else if (index_set == "cube")
    {
        return inner<Cube_iterator>(s1, s2, dirs1, dirs2);
    }
    else
    {
        std::cerr << "Invalid index set!" << std::endl;
        exit(1);
    }
}

}
