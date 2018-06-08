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

#include <boost/math/special_functions/binomial.hpp>
#include "hermite/iterators.hpp"
#include "hermite/lib.hpp"
#include "hermite/types.hpp"

using namespace std;
using boost::math::binomial_coefficient;

namespace hermite
{

u_int bissect_degree(u_int dim, u_int n_polys)
{
    u_int degree_1 = 0, degree_2 = 150;

    int img_1 = (int) binomial_coefficient<double> (degree_1 + dim, dim) - (int) n_polys;
    int img_2 = (int) binomial_coefficient<double> (degree_2 + dim, dim) - (int) n_polys;

    if (img_1 > 0 || img_2 < 0)
    {
        cout << "Can't find degree, Invalid arguments!" << endl;
        exit(0);
    }

    if (img_1 == 0)
        return degree_1;

    if (img_2 == 0)
        return degree_2;

    while (true)
    {
        u_int new_degree = (degree_1 + degree_2)/2;
        int new_img = (int) binomial_coefficient<double> (new_degree + dim, dim) - (int) n_polys;

        if (new_img < 0)
        {
            degree_1 = new_degree;
            img_1 = new_img;
        }
        else if (new_img > 0)
        {
            degree_2 = new_degree;
            img_2 = new_img;
        }
        else
        {
            return new_degree;
        }
    }
}

u_int find_dim(u_int degree, u_int n_polys)
{
    u_int dim = 0;
    u_int n_dim = (u_int) binomial_coefficient<double> (degree + dim, dim);

    while (n_dim < n_polys)
    {
        dim += 1;
        n_dim = (u_int) binomial_coefficient<double> (degree + dim, dim);
    }

    if (n_dim != n_polys)
    {
        cout << "Dimension not found, exiting..." << endl;
        exit(1);
    }

    return dim;
}

bool isAligned(const ivec & m, u_int dir)
{
    for (u_int j = 0; j < m.size(); j++)
    {
        if (m[j] != 0 && j != dir)
        {
            return false;
        }
    }
    return true;
}

bool isAligned(const ivec & m, const ivec & dirs)
{
    for (u_int j = 0; j < m.size(); j++)
    {
        bool in_dirs = false;
        for(u_int k = 0; k < dirs.size(); k++)
        {
            if (j == dirs[k])
            {
                in_dirs = true;
            }
        }

        if (m[j] != 0 && !in_dirs)
        {
            return false;
        }
    }
    return true;
}

}
