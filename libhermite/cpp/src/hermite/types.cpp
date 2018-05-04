#include "types.hpp"

namespace boost
{
    c_mat contig_mat(int rows, int cols)
    {
        auto dims = extents[rows][cols];
        return multi_array<double, 2>(dims);
    }
}
