#ifndef HERMITE_TYPES_H
#define HERMITE_TYPES_H

#include <vector>
#include <string>

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

#endif
