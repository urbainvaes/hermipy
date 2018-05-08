#ifndef CONVERTERS_H
#define CONVERTERS_H

#define WRONG_DIMENSION 1
#define ROWS_NOT_CONTIGUOUS 2

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <hermite/types.hpp>

namespace hermite
{
    boost::python::numpy::ndarray mat_to_numpy(std::mat const & input);
    boost::python::numpy::ndarray cmat_to_numpy(boost::c_mat const & input);
    boost::c_mat to_bmat(boost::python::numpy::ndarray const & input);
    std::mat to_mat(boost::python::numpy::ndarray input);
}

#endif
