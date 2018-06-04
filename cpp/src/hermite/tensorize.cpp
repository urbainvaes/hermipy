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

void set_dims(const imat & dirs, u_int & dim, ivec & dims)
{
    dim = 0;
    for (u_int i = 0; i < dirs.size(); i++)
    {
        dims[i] = dirs[i].size();
        dim += dims[i];
    }

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

vec tensorize(const mat & inputs, const imat & dirs)
{
    u_int dim;
    ivec dims(dirs.size());
    set_dims(dirs, dim, dims);
    u_int degree = bissect_degree(dims[0], inputs[0].size());
    for (u_int i = 0; i < dirs.size(); i++)
    {
        check_degree(inputs[i].size(), dims[i], degree);
    }

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
    u_int degree = inputs[0].size() - 1;
    u_int n_polys = Multi_index_iterator::size(degree, dim);
    vec results(n_polys, 0.);

    Multi_index_iterator m(dim, degree);
    u_int i,j;
    for (i = 0; !m.isFull(); i++, m.increment())
    {
        results[i] = 1;
        for (j = 0; j < dim; j++)
        {
            results[i] *= inputs[j][m[j]];
        }
    }
    return results;
}

vec tensorize(const vec & input, u_int dim, u_int dir)
{
    u_int degree = input.size() - 1;
    vec cst(degree + 1, 0.); cst[0] = 1;
    mat vecs(dim, cst);
    vecs[dir] = input;
    return tensorize(vecs);
}

template <typename T>
T tensorize(const vector<mat> & inputs, const imat & dirs)
{
    u_int dim;
    ivec dims(dirs.size());
    set_dims(dirs, dim, dims);
    u_int degree = bissect_degree(dims[0], matrix::size1(inputs[0]));
    for (u_int i = 0; i < dirs.size(); i++)
    {
        check_degree(matrix::size1(inputs[i]), dims[i], degree);
    }

    u_int n_polys = Multi_index_iterator::size(degree, dim);
    T product = matrix::construct<T>(n_polys, n_polys);

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);

    u_int i,j,k;
    for (i = 0, m1.reset(); !m1.isFull(); i++, m1.increment())
    {
        for (j = 0, m2.reset(); !m2.isFull(); j++, m2.increment())
        {
            double result = 1.;
            for (k = 0; k < dirs.size(); k++)
            {
                ivec sub1 = extract(m1.get(), dirs[k]);
                ivec sub2 = extract(m2.get(), dirs[k]);
                u_int ind1 = Multi_index_iterator::index(sub1);
                u_int ind2 = Multi_index_iterator::index(sub2);
                result *= matrix::get(inputs[k], ind1, ind2);
            }
            if (result != 0)
            {
                matrix::set(product, i, j, result);
            }
        }
    }
    return product;
}

template <typename T>
T tensorize(const vector<spmat> & inputs, const imat & dirs)
{
    u_int dim;
    ivec dims(dirs.size());
    set_dims(dirs, dim, dims);
    u_int degree = bissect_degree(dims[0], matrix::size1(inputs[0]));
    for (u_int i = 0; i < dirs.size(); i++)
    {
        check_degree(matrix::size1(inputs[i]), dims[i], degree);
    }

    u_int n_polys = Multi_index_iterator::size(degree, dim);
    T product = matrix::construct<T>(n_polys, n_polys);

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);

    u_int i,j,k;
    for (i = 0, m1.reset(); !m1.isFull(); i++, m1.increment())
    {
        for (j = 0, m2.reset(); !m2.isFull(); j++, m2.increment())
        {
            double result = 1.;
            for (k = 0; k < dirs.size(); k++)
            {
                ivec sub1 = extract(m1.get(), dirs[k]);
                ivec sub2 = extract(m2.get(), dirs[k]);
                u_int ind1 = Multi_index_iterator::index(sub1);
                u_int ind2 = Multi_index_iterator::index(sub2);
                result *= matrix::get(inputs[k], ind1, ind2);
            }
            if (result != 0)
            {
                matrix::set(product, i, j, result);
            }
        }
    }
    return product;
}

template <typename T, typename M>
T tensorize(const vector<M> & inputs)
{
    u_int dim = inputs.size();
    u_int degree = matrix::size1(inputs[0]) - 1;
    u_int n_polys = Multi_index_iterator::size(degree, dim);

    T product = matrix::construct<T>(n_polys, n_polys);

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);
    u_int i,j,k;
    for (i = 0, m1.reset(); !m1.isFull(); i++, m1.increment())
    {
        for (j = 0, m2.reset(); !m2.isFull(); j++, m2.increment())
        {
            double result = 1.;
            for (k = 0; k < dim; k++)
            {
                result *= matrix::get(inputs[k], m1[k], m2[k]);
            }
            if (result != 0)
            {
                matrix::set(product, i, j, result);
            }
        }
    }
    return product;
}


template<typename T, typename M> T
tensorize(const M & input, u_int dim, u_int dir)
{
    u_int degree = matrix::size1(input) - 1;
    T eye = matrix::eye<T>(degree + 1);
    vector<T> mats(dim, eye);
    mats[dir] = matrix::convert<T,M>(input);
    return tensorize<T>(mats);
}

template mat tensorize(const std::vector<mat> & inputs);
template mat tensorize(const std::vector<spmat> & inputs);
template spmat tensorize(const std::vector<mat> & inputs);
template spmat tensorize(const std::vector<spmat> & inputs);

template mat tensorize(const std::vector<mat> & inputs, const imat & dirs);
template mat tensorize(const std::vector<spmat> & inputs, const imat & dirs);
template spmat tensorize(const std::vector<mat> & inputs, const imat & dirs);
template spmat tensorize(const std::vector<spmat> & inputs, const imat & dirs);

template mat tensorize(const mat & input, u_int dim, u_int dir);
template mat tensorize(const spmat & input, u_int dim, u_int dir);
template spmat tensorize(const mat & input, u_int dim, u_int dir);
template spmat tensorize(const spmat & input, u_int dim, u_int dir);

}
