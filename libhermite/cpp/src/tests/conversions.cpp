#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/core/ref.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "hermite/types.hpp"
#include "hermite/python/converters.hpp"
#include "tests/conversions.hpp"

using namespace std;

namespace hermite {

boost::c_mat test3(int n)
{
    boost::c_mat result = boost::contig_mat(n, n);
    for (u_int i = 0; i < result.shape()[0]; i++)
    {
        for (u_int j = 0; j < result.shape()[1]; j++)
        {
            result[i][j] = (double) i - (double) j;
        }
    }
    return result;
}

np::ndarray test2(int n)
{
    boost::c_mat result = test3(n);
    return cmat_to_numpy(result);
}

std::mat test(int n)
{
    mat result(n, vec(n, 0.));
    for (u_int i = 0; i < result.size(); i++)
    {
        for (u_int j = 0; j < result[0].size(); j++)
        {
            result[i][j] = (double) i - (double) j;
        }
    }
    return result;
}

}
