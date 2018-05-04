#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/core/ref.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "hermite/discretize.hpp"
#include "hermite/io.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/transform.hpp"
#include "hermite/types.hpp"
#include "hermite/varf.hpp"

#include "hermite/python/converters.hpp"
#include "tests/conversions.hpp"

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

    def("test", test);
    def("test2", test2);
    def("test3", test3);
    def("to_numpy", mat_to_numpy);
    def("to_numpy", cmat_to_numpy);
    def("discretize", discretize_from_string);
    def("integrate", integrate);
    def("list_cube_indices", list_cube_indices);
    def("list_multi_indices", list_multi_indices);
    def("project", project_vec);
    def("project", project_mat);
    def("tensorize", tensorize_vecs);
    def("tensorize", tensorize_mats);
    def("tensorize", tensorize_vec);
    def("tensorize", tensorize_mat);
    def("transform", transform);
    def("triple_products", triple_products_1d);
    def("varf", varf);
    def("varfd", varfd);
}

}
