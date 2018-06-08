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

#include <iostream>
#include <assert.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "hermite/matrix.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/templates.hpp"
#include "hermite/types.hpp"
#include "hermite/lib.hpp"
#include "hermite/io.hpp"

using namespace std;

namespace hermite {

// Auxiliary functions {{{
void check_dims(const imat & dirs, u_int dim)
{
    std::vector<bool> dirs_present (dim, false);
    for (u_int i = 0; i < dirs.size(); i++)
    {
        for (u_int j = 0; j < dirs[i].size(); j++)
        {
            u_int dir = dirs[i][j];
            if (dirs_present[dir])
            {
                cout << "The same direction appears twice, exiting ..." << endl;
                exit(1);
            }
            else
            {
                dirs_present[dir] = true;
            }
        }
    }

    for (u_int i = 0; i < dirs_present.size(); i++)
    {
        if (!dirs_present[i])
        {
            cout << "Missing direction, exiting ..." << endl;
            exit(1);
        }
    }
}

void check_degree(u_int size, u_int dim, u_int degree)
{
    u_int n_polys = Multi_index_iterator::size(degree, dim);
    if (n_polys != size)
    {
        cout << "Size of input does not match dimension!" << endl;
        cout << "Expected dimension: " << n_polys << endl;
        cout << "Actual dimension: " << size << endl;
        exit(1);
    }
}

// }}}
// Tensorization of vectors {{{

vec tensorize(const mat & inputs, const imat & dirs)
{
    u_int dim = 0;
    ivec dims(dirs.size());
    for (u_int i = 0; i < dirs.size(); i++)
    {
        dims[i] = dirs[i].size();
        dim += dims[i];
    }

    u_int degree = bissect_degree(dims[0], inputs[0].size());

    #ifdef DEBUG
    check_dims(dirs, dim);
    for (u_int i = 0; i < dirs.size(); i++)
        check_degree(inputs[i].size(), dims[i], degree);
    #endif

    u_int n_polys = Multi_index_iterator::size(degree, dim);
    vec results(n_polys, 0.);
    Multi_index_iterator m(dim, degree);
    for (u_int i = 0; !m.isFull(); i++, m.increment())
    {
        results[i] = 1;
        for (u_int j = 0; j < dirs.size(); j++)
        {
            ivec sub = extract(m.get(), dirs[j]);
            u_int ind = Multi_index_iterator::index(sub);
            results[i] *= inputs[j][ind];
        }
    }
    return results;
}

vec tensorize(const mat & inputs)
{
    u_int dim = inputs.size();
    imat dirs(dim, ivec(1, 0));
    for (u_int i = 0; i < dim; i++)
        dirs[i][0] = i;
    return tensorize(inputs, dirs);
}

vec tensorize(const vec & input, u_int dim, u_int dir)
{
    u_int degree = input.size() - 1;
    vec cst(degree + 1, 0.); cst[0] = 1;
    mat vecs(dim, cst);
    vecs[dir] = input;
    return tensorize(vecs);
}

// }}}
// Tensorization of matrices {{{

void tensorize_dirs(const ivec & dA, const ivec & dB,
                    ivec & result, ivec & dA_to_ind, ivec & dB_to_ind)
{
    // Calculate compounded directions
    u_int sA = dA.size(),
          sB = dB.size();

    #ifdef DEGUB
    for (u_int i = 0; i < sA - 1; i++)
        assert(dA[i + 1] > dA[i]);

    for (u_int i = 0; i < sB - 1; i++)
        assert(dB[i + 1] > dB[i]);
    #endif

    u_int dim = sA + sB,
          iA = 0,
          iB = 0;

    dA_to_ind.resize(sA),
    dB_to_ind.resize(sB);
    result.resize(dim);

    for (u_int i = 0; i < dim; i++)
    {
        if (iB == sB || (iA < sA && dA[iA] < dB[iB]) )
        {
            result[i] = dA[iA];
            dA_to_ind[iA++] = i;
            continue;
        }
        else
        {
            result[i] = dB[iB];
            dB_to_ind[iB++] = i;
            continue;
        }
    }
}

spmat tensorize(const spmat & A, const spmat & B,
                const ivec & dA, const ivec & dB)
{
    #ifdef DEBUG
    cout << "Entering tensorize(A, B, dA, dB)" << endl;
    #endif

    // Calculate compounded directions
    u_int sA = dA.size(),
          sB = dB.size();
    u_int dim = sA + sB;

    ivec d, dA_to_ind, dB_to_ind;
    tensorize_dirs(dA, dB, d, dA_to_ind, dB_to_ind);

    u_int degree = bissect_degree(sA, A.size1());

    #ifdef DEBUG
    assert(bissect_degree(sB, B.size1()) == degree);
    #endif

    u_int n_polys = Multi_index_iterator::size(degree, dim);
    spmat product = matrix::construct<spmat>(n_polys, n_polys);

    imat multi_indices_A = list_multi_indices(sA, degree);
    imat multi_indices_B = list_multi_indices(sB, degree);

    #ifdef DEBUG
    cout << "--> Starting for loop" << endl;
    #endif
    for (cit1_t iA1 = A.begin1(); iA1 != A.end1(); ++iA1)
    {
        for (cit2_t iA2 = iA1.begin(); iA2 != iA1.end(); ++iA2)
        {
            ivec m_row_A = multi_indices_A[iA2.index1()];
            ivec m_col_A = multi_indices_A[iA2.index2()];
            double elem_A = *iA2;

            ivec multi_index_row_base(dim, 0),
                 multi_index_col_base(dim, 0);

            for (u_int i = 0; i < sA; i++)
            {
                multi_index_row_base[dA_to_ind[i]] = m_row_A[i];
                multi_index_col_base[dA_to_ind[i]] = m_col_A[i];
            }

            for (cit1_t iB1 = B.begin1(); iB1 != B.end1(); ++iB1)
            {
                for (cit2_t iB2 = iB1.begin(); iB2 != iB1.end(); ++iB2)
                {
                    ivec m_row_B = multi_indices_B[iB2.index1()];
                    ivec m_col_B = multi_indices_B[iB2.index2()];
                    double elem_B = *iB2;

                    ivec multi_index_row = multi_index_row_base,
                         multi_index_col = multi_index_col_base;

                    for (u_int i = 0; i < sB; i++)
                    {
                        multi_index_row[dB_to_ind[i]] = m_row_B[i];
                        multi_index_col[dB_to_ind[i]] = m_col_B[i];
                    }

                    u_int ind_row = Multi_index_iterator::index(multi_index_row);
                    u_int ind_col = Multi_index_iterator::index(multi_index_col);

                    if (ind_row >= matrix::size1(product) || ind_col >= matrix::size2(product))
                        continue;

                    matrix::set(product, ind_row, ind_col, elem_A * elem_B);
                }
            }
        }
    }
    return product;
}

template <typename T, typename S>
T tensorize(const vector<S> & inputs, const imat & dirs)
{
    #ifdef DEBUG
    cout << "Entering tensorize with dirs = " << endl << dirs << endl;
    cout << "--> Calculating dimensions." << endl;
    assert (inputs.size() == dirs.size());
    #endif
    u_int dim = 0;
    ivec dims(dirs.size());
    for (u_int i = 0; i < dirs.size(); i++)
    {
        dims[i] = dirs[i].size();
        dim += dims[i];
    }

    #ifdef DEBUG
    cout << "--> Dimension = " << dim << endl;
    cout << "--> Converting arguments to sparse matrices." << endl;
    #endif
    vector<spmat> sp_inputs(inputs.size());
    for (u_int i = 0; i < inputs.size(); i++)
        sp_inputs[i] = matrix::convert<spmat>(inputs[i]);

    #ifdef DEBUG
    check_dims(dirs, dim);

    u_int degree = bissect_degree(dims[0], matrix::size1(inputs[0]));
    for (u_int i = 1; i < inputs.size(); i++)
        check_degree(matrix::size1(inputs[i]), dims[i], degree);
    #endif

    spmat result = sp_inputs[0];
    ivec d = dirs[0], d_old, _a1, _a2;
    for (u_int i = 0; i < inputs.size() - 1; i++)
    {
        #ifdef DEBUG
        cout << "--> Tensorizing directions " << d << " and " << dirs[i+1] << endl;
        #endif

        d_old = d;
        tensorize_dirs(d_old, dirs[i+1], d, _a1, _a2);

        #ifdef DEBUG
        cout << "--> Tensorizing matrices" << endl;
        #endif
        result = tensorize(result, sp_inputs[i+1], d_old, dirs[i+1]);

    }

    return matrix::convert<T>(result);
}

template <typename T, typename S>
T tensorize(const vector<S> & inputs)
{
    u_int dim = inputs.size();
    imat dirs(dim, ivec(1, 0));
    for (u_int i = 0; i < dim; i++)
        dirs[i][0] = i;
    return tensorize<T,S>(inputs, dirs);

    // u_int dim = inputs.size();
    // u_int degree = matrix::size1(inputs[0]) - 1;
    // u_int n_polys = Multi_index_iterator::size(degree, dim);

    // T product = matrix::construct<T>(n_polys, n_polys);

    // Multi_index_iterator m1(dim, degree);
    // Multi_index_iterator m2(dim, degree);
    // u_int i,j,k;
    // for (i = 0, m1.reset(); !m1.isFull(); i++, m1.increment())
    // {
    //     for (j = 0, m2.reset(); !m2.isFull(); j++, m2.increment())
    //     {
    //         double result = 1.;
    //         for (k = 0; k < dim; k++)
    //         {
    //             result *= matrix::get(inputs[k], m1[k], m2[k]);
    //         }
    //         if (result != 0)
    //         {
    //             matrix::set(product, i, j, result);
    //         }
    //     }
    // }
    // return product;
}

template<typename T, typename M>
T tensorize(const M & input, u_int dim, u_int dir)
{
    u_int degree = matrix::size1(input) - 1;
    T eye = matrix::eye<T>(degree + 1);
    vector<T> mats(dim, eye);
    mats[dir] = matrix::convert<T>(input);
    return tensorize<T>(mats);
}

template mat tensorize(const std::vector<mat> & inputs, const imat & dirs);
template mat tensorize(const std::vector<spmat> & inputs, const imat & dirs);
template spmat tensorize(const std::vector<mat> & inputs, const imat & dirs);
template spmat tensorize(const std::vector<spmat> & inputs, const imat & dirs);

template mat tensorize(const std::vector<mat> & inputs);
template mat tensorize(const std::vector<spmat> & inputs);
template spmat tensorize(const std::vector<mat> & inputs);
template spmat tensorize(const std::vector<spmat> & inputs);

template mat tensorize(const mat & input, u_int dim, u_int dir);
template mat tensorize(const spmat & input, u_int dim, u_int dir);
template spmat tensorize(const mat & input, u_int dim, u_int dir);
template spmat tensorize(const spmat & input, u_int dim, u_int dir);

// }}}

}
