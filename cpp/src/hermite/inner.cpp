#include <assert.h>
#include <iostream>

#include "hermite/inner.hpp"
#include "hermite/types.hpp"
#include "hermite/iterators.hpp"
#include "hermite/lib.hpp"

namespace hermite
{

// Precondition: Dirs in increasing order
vec inner(const vec & s1,
          const vec & s2,
          const ivec & dirs1,
          const ivec & dirs2)
{
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

    u_int degree = bissect_degree(dim1, s1.size());
    assert( degree == bissect_degree(dim2, s2.size()) );

    vec result(Multi_index_iterator::size(degree, dim_result), 0.);

    Multi_index_iterator m_result(dim_result, degree);
    Multi_index_iterator m_inner(dim_inner, degree);

    for (i = 0; !m_result.isFull(); m_result.increment(), i++)
    {
        ivec m1(dim1),
             m2(dim2);

        for (u_int k = 0; k < dim_result; k++)
        {
            if (in_dirs1[k])
                m1[free_ind[k]] = m_result[k];
            else
                m2[free_ind[k]] = m_result[k];
        }

        for (j = 0, m_inner.reset(); !m_inner.isFull(); m_inner.increment(), j++)
        {
            for (u_int k = 0; k < dim_inner; k++)
            {
                m1[inner_ind1[k]] = m_inner[k];
                m2[inner_ind2[k]] = m_inner[k];
            }

            u_int ind1 = Multi_index_iterator::index(m1),
                  ind2 = Multi_index_iterator::index(m2);

            if (ind1 >= s1.size() || ind2 >= s2.size())
                continue;

            result[i] += s1[ind1] * s2[ind2];
        }
    }
    return result;
}

}
