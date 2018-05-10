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

    void set(boost::spmat & matrix, std::u_int i, std::u_int j, double val);
    void set(std::mat & matrix, std::u_int i, std::u_int j, double val);

    double get(const boost::spmat & matrix, std::u_int i, std::u_int j);
    double get(const std::mat & matrix, std::u_int i, std::u_int j);

    std::u_int size1(const std::mat &);
    std::u_int size2(const std::mat &);
    std::u_int size1(const boost::spmat &);
    std::u_int size2(const boost::spmat &);

    template<typename T> T construct(std::u_int size1, std::u_int size2);
}

#endif
