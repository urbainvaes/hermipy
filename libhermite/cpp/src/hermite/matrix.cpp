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
    void set(boost::spmat & input, std::u_int i, std::u_int j, double val)
    {
        input(i, j) = val;
    }

    void set(std::mat & input, std::u_int i, std::u_int j, double val)
    {
        input[i][j] = val;
    }

    double get(const boost::spmat & input, std::u_int i, std::u_int j)
    {
        return input(i, j);
    }

    double get(const std::mat & input, std::u_int i, std::u_int j)
    {
        return input[i][j];
    }

    template <> std::mat construct(std::u_int size1, std::u_int size2)
    {
        return std::mat(size1, std::vec(size2, 0.));
    }

    template <> boost::spmat construct(std::u_int size1, std::u_int size2)
    {
        return boost::spmat(size1, size2);
    }

    std::u_int size1(const std::mat & input)
    {
        return input.size();
    }

    std::u_int size2(const std::mat & input)
    {
        return input[0].size();
    }

    std::u_int size1(const boost::spmat & input)
    {
        return input.size1();
    }

    std::u_int size2(const boost::spmat & input)
    {
        return input.size2();
    }
}
