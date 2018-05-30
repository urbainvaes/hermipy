#include <iostream>
#include <unordered_map>

#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#include "hermite/matrix.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/templates.hpp"
#include "hermite/types.hpp"
#include "hermite/io.hpp"

using namespace std;
using boost::math::binomial_coefficient;

namespace hermite {

vec tensorize(const mat & inputs)
{
    u_int dim = inputs.size();
    u_int degree = inputs[0].size() - 1;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
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


template <typename T, typename M>
T tensorize(const vector<M> & inputs)
{
    u_int dim = inputs.size();
    u_int degree = matrix::size1(inputs[0]) - 1;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);

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

template mat tensorize(const mat & input, u_int dim, u_int dir);
template mat tensorize(const spmat & input, u_int dim, u_int dir);
template spmat tensorize(const mat & input, u_int dim, u_int dir);
template spmat tensorize(const spmat & input, u_int dim, u_int dir);

}
