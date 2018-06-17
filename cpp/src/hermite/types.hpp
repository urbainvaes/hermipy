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

#ifndef HERMITE_TYPES_H
#define HERMITE_TYPES_H

#include <vector>
#include <string>

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>


namespace hermite {

    typedef unsigned int u_int;

    // Arrays
    typedef std::vector<double> vec;
    typedef std::vector<vec> mat;
    typedef std::vector<mat> cube;

    typedef std::vector<unsigned int> ivec;
    typedef std::vector<ivec> imat;
    typedef std::vector<imat> icube;

    // Functions
    // typedef boost::function<double(vec const &)> s_func;
    typedef double (*s_func)(double*);

    // Contiguous multi-dimensional array
    typedef boost::multi_array<double, 2> cmat;
    typedef boost::numeric::ublas::matrix<double> boost_mat;

    // Sparse matrix
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major> spmat;
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>::iterator1 it1_t;
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>::iterator2 it2_t;
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>::const_iterator1 cit1_t;
    typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major>::const_iterator2 cit2_t;
}

#endif
