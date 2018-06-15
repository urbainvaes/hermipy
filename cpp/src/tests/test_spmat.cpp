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

// Compilation: g++ -I/usr/include/python3.6m -lpython3.6m -lboost_python3 -lboost_numpy3 test.cpp

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "hermite/sparse.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/templates.hpp"
#include "hermite/types.hpp"
#include "hermite/io.hpp"

#include <iostream>

using namespace std;
using namespace hermite;

hermite::mat simple(int n)
{
    hermite::mat result(n, vec(n, 0.));
    for (u_int i = 0; i < result.size(); i++)
    {
        for (u_int j = 0; j < result[0].size(); j++)
        {
            result[i][j] = (double) i*n + (double) j + 1;
        }
    }
    return result;
}

int main()
{
    int degree = 20;
    hermite::mat initial = simple(degree + 1);
    cout << initial << endl;
    hermite::spmat sp_initial = to_spmat(initial);
    vector<hermite::mat> inputs(2, initial);
    vector<hermite::spmat> sp_inputs(2, sp_initial);
    hermite::spmat sp_tensorized = tensorize<hermite::spmat, hermite::spmat>(sp_inputs, "triangle");
    hermite::mat tensorized = tensorize<hermite::mat, hermite::mat>(inputs, "triangle");
    hermite::mat difference = hermite::full(sp_tensorized) - tensorized;
    imat m = Triangle_iterator(2, degree).getList();
    for (u_int i = 0; i < tensorized.size(); i++)
    {
        for (u_int j = 0; j < tensorized.size(); j++)
        {
            double tens = initial[m[i][0]][m[j][0]] * initial[m[i][1]][m[j][1]];
            if (tens - tensorized[i][j] > 1e-12)
            {
                return 1;
            }

            if(fabs(difference[i][j]) > 1e-12)
            {
                return 1;
            }
        }
    }
    cout << "Test passed" << endl;
}
