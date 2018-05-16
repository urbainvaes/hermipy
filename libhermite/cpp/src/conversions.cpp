#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/core/ref.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "hermite/types.hpp"
#include "hermite/python/converters.hpp"
#include "tests/conversions.hpp"

using namespace std;

namespace p = boost::python;
namespace np = boost::python::numpy;
namespace hermite {

boost::cmat test3(int n)
{
    boost::cmat result = boost::contig_mat(n, n);
    for (u_int i = 0; i < result.shape()[0]; i++)
    {
        for (u_int j = 0; j < result.shape()[1]; j++)
        {
            result[i][j] = (double) i*n + (double) j + 0.1;
        }
    }
    return result;
}

np::ndarray test2(int n)
{
    boost::cmat result = test3(n);
    np::ndarray array = cmat_to_numpy(result);
    np::ndarray array2 = array;
    std::cout << p::extract<char const *>(p::str(array)) << std::endl ;
    std::cout << p::extract<char const *>(p::str(array2)) << std::endl ;
    return array;
}

std::mat test(int n)
{
    mat result(n, vec(n, 0.));
    for (u_int i = 0; i < result.size(); i++)
    {
        for (u_int j = 0; j < result[0].size(); j++)
        {
            result[i][j] = (double) i*n + (double) j;
        }
    }
    return result;
}

np::ndarray test4(int n)
{
    mat result(n, vec(n, 0.));
    for (u_int i = 0; i < result.size(); i++)
    {
        for (u_int j = 0; j < result[0].size(); j++)
        {
            result[i][j] = (double) i*n + (double) j;
        }
    }
    np::ndarray array = mat_to_numpy(result);
    std::cout << p::extract<char const *>(p::str(array)) << std::endl ;
    return array;
}


}
