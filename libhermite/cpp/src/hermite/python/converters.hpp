#ifndef CONVERTERS_H
#define CONVERTERS_H

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <hermite/types.hpp>

namespace hermite
{
    namespace np = boost::python::numpy;
    np::ndarray mat_to_numpy(std::mat const & input);
    np::ndarray cmat_to_numpy(boost::c_mat const & input);
}

#endif
