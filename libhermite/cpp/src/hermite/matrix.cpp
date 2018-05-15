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
    // Template specialization
    template <> std::mat construct(std::u_int size1, std::u_int size2)
    {
        return std::mat(size1, std::vec(size2, 0.));
    }

    template <> boost::spmat construct(std::u_int size1, std::u_int size2)
    {
        return boost::spmat(size1, size2);
    }

    template<typename T, typename S> T convert(const S & input)
    {
        return input;
    }

    template <> std::mat convert(const boost::spmat & input)
    {
        return hermite::full(input);
    }

    template <> boost::spmat convert(const std::mat & input)
    {
        return hermite::to_spmat(input);
    }

    // Template instanciation
    template boost::spmat convert(const boost::spmat & input);
    template std::mat convert(const std::mat & input);
}
