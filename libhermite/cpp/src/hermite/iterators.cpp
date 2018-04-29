#include "hermite/types.hpp"
#include "hermite/iterators.hpp"

#include <boost/math/special_functions/binomial.hpp>

namespace hermite {

    void Multi_index_iterator::increment()
    {
        if (multi_index[dim - 1] == sum)
        {
            if (sum == upper_bound)
            {
                full = true;
            }
            else
            {
                sum += 1;
                multi_index[dim - 1] = 0;
                multi_index[0] = sum;
            }
            return;
        }
        if (multi_index[0] > 0)
        {
            multi_index[0] -= 1;
            multi_index[1] += 1;
        }
        else
        {
            u_int i = 1;
            while (multi_index[i] == 0)
            {
                i += 1;
            }
            multi_index[0] += multi_index[i] - 1;
            multi_index[i] = 0;
            multi_index[i+1] += 1;
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

    std::imat list_multi_indices(u_int dim, u_int upper_bound)
    {
        std::imat result;
        for(Multi_index_iterator m(dim, upper_bound); !m.isFull(); m.increment())
        {
            result.push_back(m.get());
        }
        return result;
    }

    std::imat list_cube_indices(const std::ivec & upper_bounds)
    {
        std::imat result;
        for(Hyper_cube_iterator m(upper_bounds); !m.isFull(); m.increment())
        {
            result.push_back(m.get());
        }
        return result;
    }
}
