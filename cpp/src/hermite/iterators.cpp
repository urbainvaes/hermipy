#include "hermite/types.hpp"
#include "hermite/iterators.hpp"

#include <boost/math/special_functions/binomial.hpp>

namespace hermite
{
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

        if (i == dim - 1)
        {
            if (sum == upper_bound)
            {
                full = true;
                return;
            }
            sum += 1;
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

    void Multi_index_iterator::reset()
    {
        full = false;
        sum = 0;
        for (unsigned int i = 0; i < dim; i++)
        {
            multi_index[i] = 0;
        }
    }

    void Hyper_cube_iterator::reset()
    {
        full = false;
        for (unsigned int i = 0; i < dim; i++)
        {
            multi_index[i] = 0;
        }
    }

    u_int Multi_index_iterator::size(u_int degree, u_int dim)
    {
        using boost::math::binomial_coefficient;
        return (u_int) binomial_coefficient<double> (degree + dim, dim);
    }

    imat list_multi_indices(u_int dim, u_int upper_bound)
    {
        imat result;
        for(Multi_index_iterator m(dim, upper_bound); !m.isFull(); m.increment())
        {
            result.push_back(m.get());
        }
        return result;
    }

    imat list_cube_indices(const ivec & upper_bounds)
    {
        imat result;
        for(Hyper_cube_iterator m(upper_bounds); !m.isFull(); m.increment())
        {
            result.push_back(m.get());
        }
        return result;
    }
}
