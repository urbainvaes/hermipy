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

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost/core/ref.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "hermite/discretize.hpp"
#include "hermite/inner.hpp"
#include "hermite/iterators.hpp"
#include "hermite/project.hpp"
#include "hermite/python/converters.hpp"
#include "hermite/sparse.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/transform.hpp"
#include "hermite/types.hpp"
#include "hermite/varf.hpp"

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

    // Inner product between Hermite series
    def("inner", inner);

    // Triple products and variational formulations
    def("triple_products", triple_products_1d);

    def("varf", varf<mat>);
    def("varf_sp", varf<spmat>);
    def("varfd", varfd<mat>);
    def("varfd", varfd<spmat>);

    // Projection and tensorization of vectors
    def("project", static_cast<vec (*) (const vec & input, u_int dim, u_int dir, std::string index_set)> (& project));
    def("project", static_cast<vec (*) (const vec & input, u_int dim, const ivec & dirs, std::string index_set)> (& project));
    def("project", static_cast<mat   (*) (const mat   & input, u_int dim, u_int dir, std::string index_set)> (& project<mat>));
    def("project", static_cast<spmat (*) (const spmat & input, u_int dim, u_int dir, std::string index_set)> (& project<spmat>));
    def("project", static_cast<mat   (*) (const mat   & input, u_int dim, const ivec & dirs, std::string index_set)> (& project<mat>));
    def("project", static_cast<spmat (*) (const spmat & input, u_int dim, const ivec & dirs, std::string index_set)> (& project<spmat>));

    // Projection and tensorization of matrices
    def("tensorize", static_cast<vec (*) (const vec & input, u_int dim, u_int dir, std::string index_set)> (& tensorize));
    def("tensorize", static_cast<vec (*) (const mat & inputs, std::string index_set)> (& tensorize));
    def("tensorize", static_cast<mat (*) (const mat & inputs, u_int dim, u_int dir, std::string index_set)> (& tensorize<mat, mat>));
    def("tensorize", static_cast<mat (*) (const spmat & input, u_int dim, u_int dir, std::string index_set)> (& tensorize<mat, spmat>));
    def("tensorize", static_cast<mat (*) (const vector<mat> & input, std::string index_set)> (& tensorize<mat, mat>));
    def("tensorize", static_cast<mat (*) (const vector<spmat> & input, std::string index_set)> (& tensorize<mat, spmat>));
    def("tensorize_sp", static_cast<spmat (*) (const mat & input, u_int dim, u_int dir, std::string index_set)> (& tensorize<spmat, mat>));
    def("tensorize_sp", static_cast<spmat (*) (const spmat & input, u_int dim, u_int dir, std::string index_set)> (& tensorize<spmat, spmat>));
    def("tensorize_sp", static_cast<spmat (*) (const vector<mat> & input, std::string index_set)> (& tensorize<spmat, mat>));
    def("tensorize_sp", static_cast<spmat (*) (const vector<spmat> & input, std::string index_set)> (& tensorize<spmat, spmat>));

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
    def("list_cube_indices", Hyper_cube_iterator::list);
    def("list_multi_indices", Triangle_iterator::s_list);
}

}
