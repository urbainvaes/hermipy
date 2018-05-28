#ifndef HERMITE_TYPES_H
#define HERMITE_TYPES_H

#include <vector>
#include <string>

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace hermite {

    typedef unsigned int u_int;

    // Arrays
    typedef std::vector<double> vec;
    typedef std::vector<vec> mat;
    typedef std::vector<mat> cube;

    typedef std::vector<unsigned int> ivec;
    typedef std::vector<ivec> imat;
    typedef std::vector<imat> icube;

    // Functions
    // typedef boost::function<double(vec const &)> s_func;
    typedef double (*s_func)(double*);

    // Contiguous multi-dimensional array
    typedef boost::multi_array<double, 2> cmat;

    // Sparse matrix
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major> spmat;
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>::iterator1 it1_t;
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>::iterator2 it2_t;
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>::const_iterator1 cit1_t;
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>::const_iterator2 cit2_t;
}

#endif
