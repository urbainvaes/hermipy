#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "hermite/types.hpp"

namespace boost 
{
    c_mat contig_mat(int rows, int cols);
}

namespace matrix 
{
    inline void set(boost::spmat & input, std::u_int i, std::u_int j, double val)
    {
        input(i, j) = val;
    }

    inline void set(std::mat & input, std::u_int i, std::u_int j, double val)
    {
        input[i][j] = val;
    }

    inline std::u_int size1(const std::mat & input)
    {
        return input.size();
    }

    inline std::u_int size2(const std::mat & input)
    {
        return input[0].size();
    }

    inline std::u_int size1(const boost::spmat & input)
    {
        return input.size1();
    }

    inline std::u_int size2(const boost::spmat & input)
    {
        return input.size2();
    }

    inline double get(const boost::spmat & input, std::u_int i, std::u_int j)
    {
        return input(i, j);
    }

    inline double get(const std::mat & input, std::u_int i, std::u_int j)
    {
        return input[i][j];
    }

    template<typename T> T construct(std::u_int size1, std::u_int size2);
}

#endif
