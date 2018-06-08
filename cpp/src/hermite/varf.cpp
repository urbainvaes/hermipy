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

#include <iostream>
#include <ctime>

#include "hermite/hermite.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/transform.hpp"
#include "hermite/types.hpp"
#include "hermite/templates.hpp"
#include "hermite/varf.hpp"
#include "hermite/matrix.hpp"

#ifdef DEBUG
#include "hermite/io.hpp"
#endif

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#define MIN(i,j) (i < j ? i : j)
#define MAX(i,j) (i < j ? j : i)

using namespace std;
using boost::math::binomial_coefficient;

namespace hermite {

// Exact, does not rely on quadrature
cube triple_products_1d(int degree)
{
    cube products(2*degree + 1, mat(degree + 1, vec(degree + 1, 0.)));

    // Compute entries for i ≤ j ≤ k, and the rest by symmetry.
    int i,j,k;

    for (i = 0; i <= degree; i++)
    {
        products[i][0][i] = 1.;
    }

    vec a(2*degree + 1, 0.), b(2*degree + 1, 0.);
    for (i = 0; i <= 2*degree; i++)
    {
        a[i] = REC_A(i);
        b[i] = REC_B(i);
    }

    // i ≥ 2; j ≥ 2: h_{i+1} h_j = x REC_A(i) h_i h_j - REB_B(i) h_(i-1) h_(j)
    //                           = (REC_A(i)/REC_A(j) (h_i h_(j+1) + REC_B(j) h_i h_(j-1)) - ...
    for (i = 1; i <= degree; i++)
    {
        double ai = a[i-1];
        double bi = b[i-1];

        for (j = i; j <= degree; j++)
        {
            for (k = MAX(j - i, 0); k <= j + i - 2; k++)
                products[k][i][j] += ai/a[k] * products[k+1][i-1][j];

            for (k = MAX(j - i + 2, 0); k <= j + i - 2; k++)
                products[k][i][j] += - bi * products[k][i-2][j];

            for (k = MAX(j - i + 2, 1); k <= i + j; k++)
                products[k][i][j] += ai/a[k]*b[k] * products[k-1][i-1][j];
        }
    }

    for (i = 0; i <= degree; i++)
        for (j = i; j <= degree; j++)
            for (k = j - i; k <= j + i; k++)
                products[k][j][i] = products[k][i][j];

    return products;
}

template<typename T> T varf(
        u_int degree,
        vec const & input,
        mat const & nodes,
        mat const & weights)
{

    u_int dim = nodes.size();
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);

    #ifdef DEBUG
    cout << "Entering varf in dimension " << dim << "." << endl;
    cout << "--> Calculating Hermite transform." << endl;
    #endif

    // Hermite transform of input function
    vec Hf = transform(2*degree, input, nodes, weights, true);

    #ifdef DEBUG
    cout << "--> Determining whether to use sparse matrices." << endl;
    #endif

    // Polynomial of highest degree
    Multi_index_iterator m(dim, 2*degree); m.reset();

    #ifdef DEBUG
    u_int max_degree = 0;
    for (u_int i = 0; i < Hf.size(); i++, m.increment())
    {
        if (abs(Hf[i]) < 1e-12)
        {
            continue;
        }

        u_int sum = 0;
        for (u_int i : m.get())
        {
            sum += i;
        }

        if (sum > max_degree)
        {
            max_degree = sum;
        }
    }

    cout << "--> Maximal degree " << max_degree << endl;
    #endif

    cube products = triple_products_1d(degree);

    // To store results
    m.reset();
    T result = matrix::construct<T>(n_polys, n_polys);
    for (u_int i = 0; i < Hf.size(); i++, m.increment())
    {
        if (abs(Hf[i]) < 1e-12)
        {
            continue;
        }

        #ifdef DEBUG
        cout << "--> i = " << i << ", and m = " << m.get() << ", and Hf[i] = " << Hf[i] << endl;
        #endif

        cube factors(dim);
        for (u_int d = 0; d < dim; ++d)
        {
            factors[d] = products[m[d]];
        }

        #ifdef DEBUG
        cout << "--> Tensorizing for current mult-index." << endl;
        #endif
        auto result_iteration = tensorize<T, mat>(factors);

        #ifdef DEBUG
        cout << "--> Adding to global matrix."<< endl;
        #endif
        result = result + result_iteration*Hf[i];
    }

    #ifdef DEBUG
    cout << "--> End of varf." << endl;
    #endif

    return result;
}

template mat varf( u_int degree, vec const & input, mat const & nodes, mat const & weights);
template spmat varf( u_int degree, vec const & input, mat const & nodes, mat const & weights);

mat varfd(u_int dim, u_int degree, u_int direction, const mat & var)
{
    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);

    u_int i,j;

    mat results = mat(var.size(), vec(var.size(), 0));
    for (j = 0, m2.reset(); j < var.size(); j++, m2.increment())
    {
        if (m2[direction] == 0)
        {
           continue;
        }

        ivec diff_m2 = m2.get();
        diff_m2[direction] -= 1;
        u_int id = Multi_index_iterator::index(diff_m2);

        for (i = 0, m1.reset(); i < var.size(); i++, m1.increment())
        {
            results[i][j] = var[i][id]*sqrt(m2[direction]);
            // Entry i,j correspond to < A h_j, h_i >
        }
    }

    return results;
}

spmat varfd(u_int dim, u_int degree, u_int direction, const spmat & var)
{
    #ifdef DEBUG
    cout << "Entering varfd with sparse matrix" << endl;
    #endif

    Multi_index_iterator m(dim, degree);

    u_int i;
    imat multi_indices;
    for (i = 0, m.reset(); i < var.size1(); i++, m.increment())
    {
        multi_indices.push_back(m.get());
    }

    spmat results = spmat(var.size1(), var.size2());
    for (cit1_t i1 = var.begin1(); i1 != var.end1(); ++i1)
    {
        for (cit2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
        {
            u_int row = i2.index1();
            u_int col = i2.index2();
            ivec m_col = multi_indices[col];
            double value = *i2;

            u_int sum = 0;
            for(i = 0; i < m_col.size(); i++)
            {
                sum += m_col[i];
            }
            if (sum == degree)
            {
               continue;
            }

            ivec int_m2 = m_col;
            int_m2[direction] += 1;
            u_int id = Multi_index_iterator::index(int_m2);
            results(row, id) = value*sqrt(int_m2[direction]);
        }
    }
    return results;
}

}
