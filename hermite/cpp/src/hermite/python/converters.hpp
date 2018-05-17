#ifndef CONVERTERS_H
#define CONVERTERS_H

#define WRONG_DIMENSION 1
#define ROWS_NOT_CONTIGUOUS 2

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <hermite/types.hpp>

namespace hermite
{
    boost::python::numpy::ndarray to_numpy(const mat & input);
    boost::python::numpy::ndarray to_numpy(const cmat & input);
    cmat to_bmat(const boost::python::numpy::ndarray & input);
    mat to_mat(const boost::python::numpy::ndarray & input);
}

#endif
