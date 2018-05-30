#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/core/ref.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "hermite/discretize.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/project.hpp"
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

    class_<cmat>("contiguous_mat");

    class_<spmat>("sparse_matrix")
        .def("size1", &spmat::size1)
        .def("size2", &spmat::size2)
        ;


    // Discretization
    def("discretize", discretize_from_string);

    // Integration and Hermite transform
    def("integrate", integrate);
    def("transform", transform);

    // Triple products and variational formulations
    def("triple_products", triple_products_1d);
    def("varf", varf<mat>);
    def("varf_sp", varf<spmat>);
    def("varfd", static_cast<mat (*) (u_int dim, u_int degree, u_int dir, const mat & var)> (& varfd));
    def("varfd", static_cast<spmat (*) (u_int dim, u_int degree, u_int dir, const spmat & var)> (& varfd));


    // Projection and tensorization of vectors
    def("project", static_cast<vec (*) (const vec & input, u_int dim, u_int dir)> (& project));
    def("project", static_cast<mat   (*) (const mat   & input, u_int dim, u_int dir)> (& project<mat>));
    def("project", static_cast<spmat (*) (const spmat & input, u_int dim, u_int dir)> (& project<spmat>));
    def("project", static_cast<mat   (*) (const mat   & input, u_int dim, const ivec & dirs)> (& project<mat>));
    def("project", static_cast<spmat (*) (const spmat & input, u_int dim, const ivec & dirs)> (& project<spmat>));

    // Projection and tensorization of matrices
    def("tensorize", static_cast<vec (*) (const vec & input, u_int dim, u_int dir)> (& tensorize));
    def("tensorize", static_cast<vec (*) (const mat &)> (& tensorize));
    def("tensorize", static_cast<mat (*) (const mat & input, u_int dim, u_int dir)> (& tensorize<mat, mat>));
    def("tensorize", static_cast<mat (*) (const spmat & input, u_int dim, u_int dir)> (& tensorize<mat, spmat>));
    def("tensorize", static_cast<mat (*) (const vector<mat> & input)> (& tensorize<mat, mat>));
    def("tensorize", static_cast<mat (*) (const vector<spmat> & input)> (& tensorize<mat, spmat>));
    def("tensorize_sp", static_cast<spmat (*) (const mat & input, u_int dim, u_int dir)> (& tensorize<spmat, mat>));
    def("tensorize_sp", static_cast<spmat (*) (const spmat & input, u_int dim, u_int dir)> (& tensorize<spmat, spmat>));
    def("tensorize_sp", static_cast<spmat (*) (const vector<mat> & input)> (& tensorize<spmat, mat>));
    def("tensorize_sp", static_cast<spmat (*) (const vector<spmat> & input)> (& tensorize<spmat, spmat>));

    // Converters between data types
    def("to_numpy", static_cast<np::ndarray (*) (const mat & input)> (& to_numpy));
    def("to_numpy", static_cast<np::ndarray (*) (const cmat & input)> (& to_numpy));
    def("to_numpy", static_cast<np::ndarray (*) (const cmat & input)> (& to_numpy));
    def("to_mat", static_cast<mat (*) (const np::ndarray & input)> (& to_mat));
    def("to_bmat", to_bmat);

    // For sparse matrices
    def("row_col_val", row_col_val);
    def("to_spmat", static_cast<spmat (*) (const mat & input)> (& to_spmat));
    def("to_spmat", static_cast<spmat (*) (const vec & data, const ivec & indices, const ivec & indptr, u_int size1, u_int size2)> (& to_spmat));
    def("full", full);

    // Misc functions
    def("list_cube_indices", list_cube_indices);
    def("list_multi_indices", list_multi_indices);
}

}
