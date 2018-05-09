#ifndef HERMITE_TYPES_H
#define HERMITE_TYPES_H

#include <vector>
#include <string>

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace std {

    typedef unsigned int u_int;

    // Arrays
    typedef vector<double> vec;
    typedef vector<vec> mat;
    typedef vector<mat> cube;

    typedef vector<unsigned int> ivec;
    typedef vector<ivec> imat;
    typedef vector<imat> icube;

    // Functions
    // typedef boost::function<double(vec const &)> s_func;
    typedef double (*s_func)(double*);
}


namespace boost {

    // Contiguous multi-dimensional array
    typedef multi_array<double, 2> c_mat;
    c_mat contig_mat(int rows, int cols);
}

#endif
