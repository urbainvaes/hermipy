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

np::ndarray mat_to_numpy(mat const & input)
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
        // for (u_int j = 0; j < n_cols; j++)
        // {
        //     converted[i][j] = input[i][j];
        // }
        shape = p::make_tuple(n_cols);
        converted[i] = np::from_data(input[i].data(), dtype, shape, stride, own);
    }
    return converted;
}

np::ndarray cmat_to_numpy(boost::c_mat const & input)
{
    u_int n_rows = input.shape()[0];
    u_int n_cols = input.shape()[1];
    p::tuple shape = p::make_tuple(n_rows, n_cols);
    p::tuple strides = p::make_tuple(input.strides()[0]*sizeof(double),
                                     input.strides()[1]*sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    p::object own;
    np::ndarray converted = np::from_data(input.data(), dtype, shape, strides, own);
    return converted;
}

mat to_cpp(np::ndarray input)
{
    int dim = input.get_nd();
    Py_intptr_t const* shape = input.get_shape();
    Py_intptr_t const* strides = input.get_strides();
    double * data = reinterpret_cast<double*>(input.get_data());

    u_int n_rows = (u_int) shape[0];
    u_int n_cols = (u_int) shape[1];

    u_int stride_rows = (u_int) strides[0] / sizeof(double);
    u_int stride_cols = (u_int) strides[1] / sizeof(double);

    u_int i, j;
    mat result(n_rows);
    for (i = 0; i < n_rows; i++)
    {
        if (stride_cols == 1)
        {
            result[i] = vec(data + i*stride_rows, data + (i+1)*stride_rows);
        }
        else
        {
            cout << "Rows of the matrix must be stored contiguously" << endl; 
            throw ROWS_NOT_CONTIGUOUS;
        }
    }

    return result;

    // u_int n_rows = input.size();
    // u_int n_cols = input[0].size();
    // p::tuple shape = p::make_tuple(n_rows, n_cols);
    // p::tuple stride = p::make_tuple(sizeof(double));
    // np::dtype dtype = np::dtype::get_builtin<double>();
    // p::object own;
    // np::ndarray converted = np::zeros(shape, dtype);
}

}
