#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/core/ref.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "hermite/discretize.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/transform.hpp"
#include "hermite/types.hpp"
#include "hermite/varf.hpp"
#include "hermite/python/converters.hpp"

using namespace std;
namespace hermite {

BOOST_PYTHON_MODULE(hermite_cpp)
{
    using namespace boost::python;

    // Initialize numpy
    Py_Initialize();
    boost::python::numpy::initialize();

    class_<vec>("double_vec")
        .def(vector_indexing_suite<vec>())
        ;

    class_<mat>("double_mat")
        .def(vector_indexing_suite<mat>())
        ;

    class_<cube>("double_cube")
        .def(vector_indexing_suite<cube>())
        ;

    class_<ivec>("int_vec")
        .def(vector_indexing_suite<ivec>())
        ;

    class_<imat>("int_mat")
        .def(vector_indexing_suite<imat>())
        ;

    class_<boost::c_mat>("contiguous_mat");


    // Core functions
    def("discretize", discretize_from_string);
    def("integrate", integrate);
    def("transform", transform);
    def("triple_products", triple_products_1d);
    def("varf", varf);
    def("varfd", varfd);

    // Tensorization
    def("project", static_cast<std::vec (*) (const std::vec & input, std::u_int dim, std::u_int dir)> (& project));
    def("project", static_cast<std::mat (*) (const std::mat & input, std::u_int dim, std::u_int dir)> (& project));
    def("tensorize", static_cast<std::vec (*) (const std::mat &)> (& tensorize));
    def("tensorize", static_cast<std::mat (*) (const std::cube &)> (& tensorize));
    def("tensorize", static_cast<std::vec (*) (const std::vec & input, std::u_int dim, std::u_int dir)> (& tensorize));
    def("tensorize", static_cast<std::mat (*) (const std::mat & input, std::u_int dim, std::u_int dir)> (& tensorize));
    def("tensorize", static_cast<boost::spmat (*) (const std::vector<boost::spmat> &)>(& tensorize));

    // Converters between data types
    def("to_numpy", mat_to_numpy);
    def("to_numpy", cmat_to_numpy);
    def("to_mat", to_mat);
    def("to_bmat", to_bmat);

    // Misc functions
    def("list_cube_indices", list_cube_indices);
    def("list_multi_indices", list_multi_indices);
}

}
