#ifndef CONVERTERS_H
#define CONVERTERS_H

#define WRONG_DIMENSION 1
#define ROWS_NOT_CONTIGUOUS 2

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <hermite/types.hpp>

namespace hermite
{
    boost::python::numpy::ndarray to_numpy(const std::mat & input);
    boost::python::numpy::ndarray to_numpy(const boost::c_mat & input);
    boost::c_mat to_bmat(const boost::python::numpy::ndarray & input);
    std::mat to_mat(const boost::python::numpy::ndarray & input);
}

#endif
