#include "hermite/matrix.hpp"

namespace boost
{
    c_mat contig_mat(int rows, int cols)
    {
        auto dims = extents[rows][cols];
        return multi_array<double, 2>(dims);
    }
}

namespace matrix
{
    template <> std::mat construct(std::u_int size1, std::u_int size2)
    {
        return std::mat(size1, std::vec(size2, 0.));
    }

    template <> boost::spmat construct(std::u_int size1, std::u_int size2)
    {
        return boost::spmat(size1, size2);
    }
}
