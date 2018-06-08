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

#ifndef CONVERTERS_H
#define CONVERTERS_H

#define WRONG_DIMENSION 1
#define ROWS_NOT_CONTIGUOUS 2

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <hermite/types.hpp>

namespace hermite
{
    boost::python::numpy::ndarray to_numpy(const mat & input);
    boost::python::numpy::ndarray to_numpy(const cmat & input);
    cmat to_bmat(const boost::python::numpy::ndarray & input);
    mat to_mat(const boost::python::numpy::ndarray & input);
}

#endif
