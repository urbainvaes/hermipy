// Compilation: g++ -I/usr/include/python3.6m -lpython3.6m -lboost_python3 -lboost_numpy3 test.cpp

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "boost/multi_array.hpp"

#include <iostream>

using namespace std;

namespace p = boost::python;
namespace np = boost::python::numpy;

np::ndarray test_array()
{
    double vector[] = {0., 1., 2., 3., 4., 5.};
    p::tuple shape = p::make_tuple(6);
    np::dtype dtype = np::dtype::get_builtin<double>();
    p::tuple stride = p::make_tuple(sizeof(double));
    np::ndarray result = np::from_data(vector, dtype, shape, stride, p::object());
    cout << p::extract<char const *>(p::str(result)) << endl;
    return result.copy();
}

typedef boost::multi_array<double, 2> mat;
np::ndarray test_multi_array(int n)
{
    double value = 0.;
    auto dims = boost::extents[n][n];
    mat array = boost::multi_array<double, 2>(dims);
    for (unsigned int i = 0; i < array.shape()[0]; i++)
        for (unsigned int j = 0; j < array.shape()[1]; j++)
            array[i][j] = ++value;
    p::tuple shape = p::make_tuple(n, n);
    p::tuple strides = p::make_tuple(array.strides()[0]*sizeof(double),
                                     array.strides()[1]*sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    np::ndarray result = np::from_data(array.data(), dtype, shape, strides, p::object());
    cout << p::extract<char const *>(p::str(result)) << endl;
    return result.copy();
}

np::ndarray convert(mat const & array)
{
    p::tuple shape = p::make_tuple(array.shape()[0], array.shape()[1]);
    p::tuple strides = p::make_tuple(array.strides()[0]*sizeof(double),
                                     array.strides()[1]*sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    np::ndarray result = np::from_data(array.data(), dtype, shape, strides, p::object());
    cout << p::extract<char const *>(p::str(result)) << endl;
    return result.copy();
}

typedef boost::multi_array<double, 2> mat;
mat test_return_multi_array(int n)
{
    double value = 0.;
    auto dims = boost::extents[n][n];
    mat array = boost::multi_array<double, 2>(dims);
    for (unsigned int i = 0; i < array.shape()[0]; i++)
        for (unsigned int j = 0; j < array.shape()[1]; j++)
            array[i][j] = ++value;
    return mat;
}

int main()
{
    Py_Initialize();
    boost::python::numpy::initialize();

    np::ndarray test1 = test_array();
    cout << p::extract<char const *>(p::str(test1)) << endl;

    int n = 5;
    auto dims = boost::extents[n][n];
    mat array = mat(dims); double value = 0.;
    for (unsigned int i = 0; i < array.shape()[0]; i++)
        for (unsigned int j = 0; j < array.shape()[1]; j++)
            array[i][j] = ++value;
    np::ndarray test2 = convert(array);
    test2 = test_multi_array(n);
    cout << p::extract<char const *>(p::str(test2)) << endl;
}
