#include "hermite/types.hpp"
#include "hermite/iterators.hpp"

using namespace hermite;

void Multi_index_iterator::increment()
{
    unsigned int i = dim - 1;
    while(sum == upper_bound && i > 0) {
        sum -= multi_index[i];
        multi_index[i] = 0;
        i -= 1;
    }
    if (i == 0 && sum == upper_bound)
    {
        full = true;
    }
    else
    {
        multi_index[i] += 1;
        sum += 1;
    }
}

void Hyper_cube_iterator::increment()
{
    unsigned int i = dim - 1;
    while(multi_index[i] == upper_bounds[i] - 1 && i > 0) {
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
