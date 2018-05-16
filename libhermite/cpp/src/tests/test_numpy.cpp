#include <cmath>
#include <iostream>

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "boost/multi_array.hpp"
#include "boost/multi_array.hpp"

#include "hermite/python/converters.hpp"
#include "hermite/types.hpp"

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
    return result.copy();
}

np::ndarray convert(hermite::cmat const & array)
{
    p::tuple shape = p::make_tuple(array.shape()[0], array.shape()[1]);
    p::tuple strides = p::make_tuple(array.strides()[0]*sizeof(double),
                                     array.strides()[1]*sizeof(double));
    np::dtype dtype = np::dtype::get_builtin<double>();
    np::ndarray result = np::from_data(array.data(), dtype, shape, strides, p::object());
    cout << p::extract<char const *>(p::str(result)) << endl;
    return result.copy();
}

hermite::cmat test_return_multi_array(int n)
{
    double value = 0.;
    auto dims = boost::extents[n][n];
    hermite::cmat array = boost::multi_array<double, 2>(dims);
    for (unsigned int i = 0; i < array.shape()[0]; i++)
        for (unsigned int j = 0; j < array.shape()[1]; j++)
            array[i][j] = ++value;
    return array;
}

int main()
{
    Py_Initialize();
    boost::python::numpy::initialize();
    u_int i,j;

    np::ndarray narray = test_array();
    cout << p::extract<char const *>(p::str(narray)) << endl;

    u_int n = 5;
    auto dims = boost::extents[n][n];
    hermite::cmat array = hermite::cmat(dims); double value = 0.;
    for (unsigned int i = 0; i < array.shape()[0]; i++)
        for (unsigned int j = 0; j < array.shape()[1]; j++)
            array[i][j] = ++value;

    hermite::cmat returned = test_return_multi_array(n);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if(fabs(returned[i][j] - array[i][j]) > 1e-12)
                return 1;

    hermite::mat converted = hermite::to_mat(hermite::to_numpy(array));
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if(fabs(converted[i][j] - array[i][j]) > 1e-12)
                return 1;

    cout << "Test passed" << endl;
}
