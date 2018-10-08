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
    class_<boost_mat>("boost_mat")
        .def("size1", &spmat::size1)
        .def("size2", &spmat::size2)
        ;

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
    def("triple_products_fourier", triple_products_fourier);

    def("varf", varf<boost_mat>);
    def("varf_sp", varf<spmat>);
    def("varfd", varfd<boost_mat>);
    def("varfd", varfd<spmat>);

    // Projection and tensorization of vectors
    def("project", & project_vec_1d);
    def("project", & project_vec_nd);
    def("project", & project_mat_1d<mat>);
    def("project", & project_mat_1d<spmat>);
    def("project", & project_mat_nd<mat>);
    def("project", & project_mat_nd<spmat>);

    // Projection and tensorization of matrices
    def("tensorize", & tensorize_vec_id);
    def("tensorize", & tensorize_vecs_axes);
    def("tensorize", & tensorize_vecs_dirs);

    def("tensorize", & tensorize_mat_id<mat,mat>);
    def("tensorize", & tensorize_mat_id<mat,spmat>);
    def("tensorize", & tensorize_mats_axes<mat,mat>);
    def("tensorize", & tensorize_mats_axes<mat,spmat>);
    def("tensorize", & tensorize_mats_dirs<mat,mat>);
    def("tensorize", & tensorize_mats_dirs<mat,spmat>);

    def("tensorize_sp", & tensorize_mat_id<spmat,mat>);
    def("tensorize_sp", & tensorize_mat_id<spmat,spmat>);
    def("tensorize_sp", & tensorize_mats_axes<spmat,mat>);
    def("tensorize_sp", & tensorize_mats_axes<spmat,spmat>);
    def("tensorize_sp", & tensorize_mats_dirs<spmat,mat>);
    def("tensorize_sp", & tensorize_mats_dirs<spmat,spmat>);

    // Converters between data types
    def("to_numpy", & to_numpy<mat>);
    def("to_numpy", & to_numpy<cmat>);
    def("to_numpy", & to_numpy<boost_mat>);

    def("to_mat", static_cast<mat (*) (const np::ndarray & input)> (& to_mat));
    def("to_bmat", to_bmat);
    def("to_boost_mat", to_boost_mat);

    // For sparse matrices
    def("row_col_val", row_col_val);
    def("to_spmat", static_cast<spmat (*) (const mat & input)> (& to_spmat));
    def("to_spmat", static_cast<spmat (*) (const vec & data, const ivec & indices, const ivec & indptr, u_int size1, u_int size2)> (& to_spmat));
    def("full", full);

    // Iterator functions
    def("list_cube_indices", Grid_iterator::list);

    def("triangle_index", Triangle_iterator::s_index);
    def("cube_index", Cube_iterator::s_index);

    def("triangle_list_indices", Triangle_iterator::s_list);
    def("cross_list_indices", Cross_iterator::s_list);
    def("cross_nc_list_indices", Cross_iterator_nc::s_list);
    def("cube_list_indices", Cube_iterator::s_list);
    def("rectangle_list_indices", Rectangle_iterator::s_list);

    def("triangle_get_degree", Triangle_iterator::s_get_degree);
    def("cross_get_degree", Cross_iterator::s_get_degree);
    def("cross_nc_get_degree", Cross_iterator_nc::s_get_degree);
    def("cube_get_degree", Cube_iterator::s_get_degree);
    def("rectangle_get_degree", Rectangle_iterator::s_get_degree);

    def("triangle_size", Triangle_iterator::s_size);
    def("cross_size", Cross_iterator::s_size);
    def("cross_nc_size", Cross_iterator_nc::s_size);
    def("cube_size", Cube_iterator::s_size);
    def("rectangle_size", Rectangle_iterator::s_size);
}

}
