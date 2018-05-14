#include <iostream>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/core/ref.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "hermite/python/converters.hpp"

using namespace std;

namespace hermite {
namespace p = boost::python;
namespace np = boost::python::numpy;

np::ndarray to_numpy(const mat & input)
{
    u_int n_rows = input.size();
    u_int n_cols = input[0].size();
    p::tuple shape = p::make_tuple(n_rows, n_cols);
    p::tuple stride = p::make_tuple(sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    p::object own;
    np::ndarray converted = np::zeros(shape, dtype);

    for (u_int i = 0; i < n_rows; i++)
    {
        shape = p::make_tuple(n_cols);
        converted[i] = np::from_data(input[i].data(), dtype, shape, stride, own);
    }
    return converted.copy();
}

np::ndarray to_numpy(const boost::c_mat & input)
{
    u_int n_rows = input.shape()[0];
    u_int n_cols = input.shape()[1];
    p::tuple shape = p::make_tuple(n_rows, n_cols);
    p::tuple strides = p::make_tuple(input.strides()[0]*sizeof(double),
                                     input.strides()[1]*sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    np::ndarray converted = np::from_data(input.data(), dtype, shape, strides, p::object());
    return converted.copy();
}


boost::c_mat to_bmat(const np::ndarray & input)
{
    if (input.get_nd() != 2)
    {
        cout << "Dimension must be 2" << endl;
        throw WRONG_DIMENSION;
    }

    Py_intptr_t const* shape = input.get_shape();
    Py_intptr_t const* strides = input.get_strides();
    double * data = reinterpret_cast<double*>(input.get_data());

    u_int n_rows = (u_int) shape[0];
    u_int n_cols = (u_int) shape[1];

    u_int stride_cols = (u_int) strides[1] / sizeof(double);

    if (stride_cols != 1)
    {
        cout << "Rows of the matrix must be stored contiguously" << endl;
        throw ROWS_NOT_CONTIGUOUS;
    }

    auto sizes = boost::extents[n_rows][n_cols];
    boost::c_mat result = boost::c_mat(boost::multi_array_ref<double, 2>(data, sizes)); // copy
    return result;
}


mat to_mat(const np::ndarray & input)
{
    if (input.get_nd() != 2)
    {
        cout << "Dimension must be 2" << endl;
        throw WRONG_DIMENSION;
    }

    Py_intptr_t const* shape = input.get_shape();
    Py_intptr_t const* strides = input.get_strides();
    double * data = reinterpret_cast<double*>(input.get_data());

    u_int n_rows = (u_int) shape[0];

    u_int stride_rows = (u_int) strides[0] / sizeof(double);
    u_int stride_cols = (u_int) strides[1] / sizeof(double);

    if (stride_cols != 1)
    {
        cout << "Rows of the matrix must be stored contiguously" << endl;
        throw ROWS_NOT_CONTIGUOUS;
    }

    mat result(n_rows);
    for (u_int i = 0; i < n_rows; i++)
    {
        result[i] = vec(data + i*stride_rows, data + (i+1)*stride_rows);
    }

    return result;
}

}
