#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "tests/conversions.hpp"
#include "hermite/python/converters.hpp"
#include "hermite/io.hpp"

#include <iostream>

using namespace std;
using namespace hermite;

namespace p = boost::python;
namespace np = boost::python::numpy;

int main()
{
    Py_Initialize();
    boost::python::numpy::initialize();

    int n = 1000;
    boost::c_mat result = test3(n);
    np::ndarray array_cmat = cmat_to_numpy(result);
    cout << array_cmat.strides(0) << endl;

    // mat result_mat = test(n);
    // // hermite::printMat(result_mat, " ");
    // np::ndarray array = hermite::mat_to_numpy(result_mat);


    // mat reconv = hermite::to_cpp(array);
    // cout << p::extract<char const *>(p::str(array_cmat)) << endl;
    // cout << p::extract<char const *>(p::str(array)) << endl;

    // hermite::printMat(reconv, " ");
    // cout << "Hello, world!" << endl;
}
