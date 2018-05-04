#ifndef CONVERSIONS_H
#define CONVERSIONS_H

#include "hermite/types.hpp"

namespace p = boost::python;
namespace np = boost::python::numpy;

namespace hermite {

    boost::c_mat test3(int n);
    np::ndarray test2(int n);
    std::mat test(int n);

}

#endif
