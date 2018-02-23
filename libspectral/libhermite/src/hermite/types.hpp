#ifndef HERMITE_TYPES_H
#define HERMITE_TYPES_H

#include <vector>

namespace hermite {
    // Arrays
    typedef std::vector<double> vec;
    typedef std::vector<vec> mat;
    typedef std::vector<mat> cube;

    typedef std::vector<unsigned int> ivec;

    // Functions
    // typedef boost::function<double(vec const &)> s_func;
    typedef double (*s_func)(double*);
}

#endif
