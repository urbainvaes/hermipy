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
#include "hermite/sparse.hpp"
#include "hermite/python/converters.hpp"

using namespace std;

namespace p = boost::python;
namespace np = boost::python::numpy;
namespace bnu = boost::numeric::ublas;

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

    class_<boost::spmat>("sparse_matrix")
        .def("size1", &boost::spmat::size1)
        .def("size2", &boost::spmat::size2)
        ;


    // Discretization
    def("discretize", discretize_from_string);

    // Integration and Hermite transform
    def("integrate", integrate);
    def("transform", transform);

    // Triple products and variational formulations
    def("triple_products", triple_products_1d);
    def("varf", varf<std::mat>);
    def("varfd", varfd);

    // Sparse versions
    def("varf_sp", varf<boost::spmat>);

    // Projection and tensorization of vectors
    def("project", static_cast<std::vec (*) (const std::vec & input, std::u_int dim, std::u_int dir)> (& project));
    def("tensorize", static_cast<std::vec (*) (const std::vec & input, std::u_int dim, std::u_int dir)> (& tensorize));
    def("tensorize", static_cast<std::vec (*) (const std::mat &)> (& tensorize));

    // Projection and tensorization of matrices
    def("project", static_cast<std::mat (*) (const std::mat & input, std::u_int dim, std::u_int dir)> (& project));
    def("tensorize", static_cast<std::mat (*) (const std::mat & input, std::u_int dim, std::u_int dir)> (& tensorize<std::mat, std::mat>));
    def("tensorize", static_cast<std::mat (*) (const boost::spmat & input, std::u_int dim, std::u_int dir)> (& tensorize<std::mat, boost::spmat>));
    def("tensorize", static_cast<std::mat (*) (const std::vector<std::mat> & input)> (& tensorize<std::mat, std::mat>));
    def("tensorize", static_cast<std::mat (*) (const std::vector<boost::spmat> & input)> (& tensorize<std::mat, boost::spmat>));
    def("tensorize_sp", static_cast<boost::spmat (*) (const std::mat & input, std::u_int dim, std::u_int dir)> (& tensorize<boost::spmat, std::mat>));
    def("tensorize_sp", static_cast<boost::spmat (*) (const boost::spmat & input, std::u_int dim, std::u_int dir)> (& tensorize<boost::spmat, boost::spmat>));
    def("tensorize_sp", static_cast<boost::spmat (*) (const std::vector<std::mat> & input)> (& tensorize<boost::spmat, std::mat>));
    def("tensorize_sp", static_cast<boost::spmat (*) (const std::vector<boost::spmat> & input)> (& tensorize<boost::spmat, boost::spmat>));

    // Converters between data types
    def("to_numpy", static_cast<np::ndarray (*) (const std::mat & input)> (& to_numpy));
    def("to_numpy", static_cast<np::ndarray (*) (const boost::c_mat & input)> (& to_numpy));
    def("to_numpy", static_cast<np::ndarray (*) (const boost::c_mat & input)> (& to_numpy));
    def("to_mat", static_cast<std::mat (*) (const np::ndarray & input)> (& to_mat));
    def("to_bmat", to_bmat);

    // For sparse matrices
    def("row_col_val", row_col_val);
    def("to_spmat", to_spmat);
    def("full", full);

    // Misc functions
    def("list_cube_indices", list_cube_indices);
    def("list_multi_indices", list_multi_indices);
}

}
