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

#include <cmath>
#include <functional>
#include <numeric>
#include <iostream>
#include <fstream>
#include <tuple>
#include <utility>
#include <vector>

#include <boost/function.hpp>
#include <boost/core/ref.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "hermite/types.hpp"
#include "hermite/hermite.hpp"
#include "hermite/iterators.hpp"
#include "hermite/transform.hpp"

using namespace std;

namespace hermite {

vec transform(
        u_int degree,
        vec const & input,
        mat const & nodes,
        mat const & weights,
        bool forward)
{
    u_int dim = nodes.size();
    using boost::math::binomial_coefficient;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);

    u_int i,j,k;

    u_int n_points_tot = 1;
    ivec n_points(dim);
    for (i = 0; i < dim; i++)
    {
        n_points[i] = nodes[i].size();
        n_points_tot *= n_points[i];
    }

    vec output;
    if (forward)
    {
        output = vec(n_polys, 0);
    }
    else
    {
        output = vec(n_points_tot, 0);
    }

    // Compute Hermite polynomials in each dimension
    cube herm_vals_1d(dim);
    for (i = 0; i < dim; ++i)
    {
        herm_vals_1d[i] = mat(n_points[i], vec(degree + 1, 0));
        for (j = 0; j < n_points[i]; j++)
        {
            hermite_eval(nodes[i][j], degree, herm_vals_1d[i][j]);
        }
    }

    Hyper_cube_iterator p(n_points);
    for (i = 0; i < n_points_tot; i++, p.increment())
    {
        double weight = 1;
        if (forward)
        {
            for (j = 0; j < dim; j++)
            {
                weight *= weights[j][p[j]];
            }
        }

        Multi_index_iterator m(dim, degree);
        for (j = 0; j < n_polys; j++, m.increment())
        {
            double val_at_point = 1;
            for (k = 0; k < dim; k++)
            {
                if (m[k] != 0)
                {
                    val_at_point *= herm_vals_1d[k][p[k]][m[k]];
                }
            }
            if(forward)
            {
                output[j] += input[i] * weight * val_at_point;
            }
            else
            {
                output[i] += input[j] * weight * val_at_point;
            }
        }
    }

    return output;
}

double integrate(
        vec const & f_grid,
        mat const & nodes,
        mat const & weights)
{
    vec integral = transform(0, f_grid, nodes, weights, true);
    return integral[0];
}

}
