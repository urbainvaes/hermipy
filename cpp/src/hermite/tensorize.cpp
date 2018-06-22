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

template<typename Iterator>
void check_degree(u_int size, u_int dim, u_int degree)
{
    u_int n_polys = Iterator::s_size(dim, degree);
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

template<typename Iterator>
vec tensorize(const mat & inputs, const imat & dirs)
{
    u_int dim = 0;
    ivec dims(dirs.size());
    for (u_int i = 0; i < dirs.size(); i++)
    {
        dims[i] = dirs[i].size();
        dim += dims[i];
    }

    u_int degree = Iterator::s_get_degree(dims[0], inputs[0].size());
    std::vector<Iterator> it_mul;
    for (u_int i = 0; i < dirs.size(); i++)
        it_mul.push_back(Iterator(dims[i], degree));

    #ifdef DEBUG
    check_dims(dirs, dim);
    for (u_int i = 0; i < dirs.size(); i++)
        check_degree<Iterator>(inputs[i].size(), dims[i], degree);
    #endif

    u_int n_polys = Iterator::s_size(dim, degree);
    vec results(n_polys, 0.);
    Iterator m(dim, degree);
    for (u_int i = 0; !m.isFull(); i++, m.increment())
    {
        results[i] = 1;
        for (u_int j = 0; j < dirs.size(); j++)
        {
            ivec sub = extract(m.get(), dirs[j]);
            u_int ind = it_mul[j].index(sub);
            results[i] *= inputs[j][ind];
        }
    }
    return results;
}

vec tensorize_vecs_dirs(const mat & inputs, const imat & dirs, std::string index_set)
{
    if (index_set == "cross")
    {
        return tensorize<Cross_iterator>(inputs, dirs);
    }
    else if (index_set == "triangle")
    {
        return tensorize<Triangle_iterator>(inputs, dirs);
    }
    else if (index_set == "cube")
    {
        return tensorize<Cube_iterator>(inputs, dirs);
    }
    else
    {
        std::cerr << "Invalid index set!" << std::endl;
        exit(1);
    }
}

vec tensorize_vecs_axes(const mat & inputs, std::string index_set)
{
    u_int dim = inputs.size();
    imat dirs(dim, ivec(1, 0));
    for (u_int i = 0; i < dim; i++)
        dirs[i][0] = i;
    return tensorize_vecs_dirs(inputs, dirs, index_set);
}

vec tensorize_vec_id(const vec & input, u_int dim, u_int dir, std::string index_set)
{
    u_int degree = input.size() - 1;
    vec cst(degree + 1, 0.); cst[0] = 1;
    mat vecs(dim, cst);
    vecs[dir] = input;
    return tensorize_vecs_axes(vecs, index_set);
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

template<typename Iterator>
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

    u_int degree = Iterator::s_get_degree(sA, A.size1());

    #ifdef DEBUG
    assert(Iterator::s_get_degree(sB, B.size1()) == degree);
    #endif

    u_int n_polys = Iterator::s_size(dim, degree);
    spmat product = matrix::construct<spmat>(n_polys, n_polys);

    Iterator it_product = Iterator(dim, degree);
    imat multi_indices_A = Iterator::s_list(sA, degree),
         multi_indices_B = Iterator::s_list(sB, degree);

    #ifdef DEBUG
    cout << "--> With the following multi-index sets:" << endl;
    cout << multi_indices_A << endl << multi_indices_B << endl;;
    #endif

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

                    bool has_row = it_product.has(multi_index_row),
                         has_col = it_product.has(multi_index_col);

                    if (!has_row || !has_col)
                        continue;

                    u_int ind_row = it_product.index(multi_index_row),
                          ind_col = it_product.index(multi_index_col);

                    matrix::set(product, ind_row, ind_col, elem_A * elem_B);
                }
            }
        }
    }
    return product;
}

template <typename Iterator, typename T, typename S>
T _tensorize_mats_dirs(const vector<S> & inputs, const imat & dirs)
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

    u_int degree = Iterator::s_get_degree(dims[0], matrix::size1(inputs[0]));
    for (u_int i = 1; i < inputs.size(); i++)
        check_degree<Iterator>(matrix::size1(inputs[i]), dims[i], degree);
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
        result = tensorize<Iterator>(result, sp_inputs[i+1], d_old, dirs[i+1]);
    }

    return matrix::convert<T>(result);
}

template <typename Iterator, typename T, typename S>
T _tensorize_mats_axes(const vector<S> & inputs)
{
    u_int dim = inputs.size();
    imat dirs(dim, ivec(1, 0));
    for (u_int i = 0; i < dim; i++)
        dirs[i][0] = i;
    return _tensorize_mats_dirs<Iterator,T,S>(inputs, dirs);
}

template <typename T, typename S>
T tensorize_mats_dirs(const vector<S> & inputs, const imat & dirs, std::string index_set)
{
    if (index_set == "cross")
    {
        return _tensorize_mats_dirs<Cross_iterator,T,S>(inputs, dirs);
    }
    else if (index_set == "triangle")
    {
        return _tensorize_mats_dirs<Triangle_iterator,T,S>(inputs, dirs);
    }
    else if (index_set == "cube")
    {
        return _tensorize_mats_dirs<Cube_iterator,T,S>(inputs, dirs);
    }
    else
    {
        std::cerr << "Invalid index set!" << std::endl;
        exit(1);
    }
}

template <typename T, typename S>
T tensorize_mats_axes(const vector<S> & inputs, std::string index_set)
{
    if (index_set == "cross")
    {
        return _tensorize_mats_axes<Cross_iterator,T,S>(inputs);
    }
    else if (index_set == "triangle")
    {
        return _tensorize_mats_axes<Triangle_iterator,T,S>(inputs);
    }
    else if (index_set == "cube")
    {
        return _tensorize_mats_axes<Cube_iterator,T,S>(inputs);
    }
    else
    {
        std::cerr << "Invalid index set!" << std::endl;
        exit(1);
    }
}

template<typename T, typename S>
T tensorize_mat_id(const S & input, u_int dim, u_int dir, std::string index_set)
{
    u_int degree = matrix::size1(input) - 1;
    T eye = matrix::eye<T>(degree + 1);
    vector<T> mats(dim, eye);
    mats[dir] = matrix::convert<T>(input);
    return tensorize_mats_axes<T>(mats, index_set);
}

template mat _tensorize_mats_axes<Cube_iterator>(const std::vector<mat> & input);
template mat _tensorize_mats_axes<Cross_iterator>(const std::vector<mat> & input);
template mat _tensorize_mats_axes<Triangle_iterator>(const std::vector<mat> & input);
template spmat _tensorize_mats_axes<Cube_iterator>(const std::vector<mat> & input);
template spmat _tensorize_mats_axes<Cross_iterator>(const std::vector<mat> & input);
template spmat _tensorize_mats_axes<Triangle_iterator>(const std::vector<mat> & input);
template boost_mat _tensorize_mats_axes<Cube_iterator>(const std::vector<mat> & input);
template boost_mat _tensorize_mats_axes<Cross_iterator>(const std::vector<mat> & input);
template boost_mat _tensorize_mats_axes<Triangle_iterator>(const std::vector<mat> & input);

template mat tensorize_mats_dirs(const std::vector<mat> & inputs, const imat & dir, std::string index_sets);
template mat tensorize_mats_dirs(const std::vector<spmat> & inputs, const imat & dir, std::string index_sets);
template spmat tensorize_mats_dirs(const std::vector<mat> & inputs, const imat & dir, std::string index_sets);
template spmat tensorize_mats_dirs(const std::vector<spmat> & inputs, const imat & dir, std::string index_sets);

template mat tensorize_mats_axes(const std::vector<mat> & input, std::string index_sets);
template mat tensorize_mats_axes(const std::vector<spmat> & input, std::string index_sets);
template spmat tensorize_mats_axes(const std::vector<mat> & input, std::string index_sets);
template spmat tensorize_mats_axes(const std::vector<spmat> & input, std::string index_sets);

template mat tensorize_mat_id(const mat & input, u_int dim, u_int dir, std::string index_set);
template mat tensorize_mat_id(const spmat & input, u_int dim, u_int dir, std::string index_set);
template spmat tensorize_mat_id(const mat & input, u_int dim, u_int dir, std::string index_set);
template spmat tensorize_mat_id(const spmat & input, u_int dim, u_int dir, std::string index_set);

// }}}

}
