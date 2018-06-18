/*
 * Copyright (C) 2018 Urbain Vaes
 *
 * This file is part of hermipy, a python/C++ library for automating the
 * Hermite Galerkin method.
 *
 * hermipy is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hermipy is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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

np::ndarray to_numpy(const cmat & input)
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


cmat to_bmat(const np::ndarray & input)
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
    cmat result = cmat(boost::multi_array_ref<double, 2>(data, sizes)); // copy
    return result;
}

boost_mat to_boost_mat(const np::ndarray & input)
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

    boost_mat result (n_rows, n_cols);
    for (u_int i = 0; i < n_rows; ++i)
        for (u_int j = 0; j < n_cols; j++)
            result(i, j) = data[i*n_rows + j];

    return result;
}

np::ndarray boost_to_numpy(const boost_mat & input)
{
    u_int n_rows = input.size1();
    u_int n_cols = input.size2();
    p::tuple shape = p::make_tuple(n_rows, n_cols);
    p::tuple strides = p::make_tuple(n_cols * sizeof(double), sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    np::ndarray converted = np::from_data(input.data().begin(), dtype, shape, strides, p::object());
    #ifdef DEBUG
    std::cout << "Converted matrix:" << std::endl;
    std::cout << p::extract<char const *>(p::str(converted)) << std::endl;
    #endif
    return converted.copy();
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
