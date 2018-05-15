#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <iostream>

#include "hermite/matrix.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/templates.hpp"
#include "hermite/types.hpp"

using namespace std;
using boost::math::binomial_coefficient;

namespace hermite {

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

vec project(const vec & input, u_int dim, u_int dir)
{
    u_int degree = 0, n_polys = input.size();
    while ((u_int) binomial_coefficient<double> (degree + dim, dim) != n_polys)
    {
        degree++;
    }
    vec results(degree + 1, 0.);

    Multi_index_iterator m(dim, degree);
    for (u_int i = 0; !m.isFull(); i++, m.increment())
    {
        if (isAligned(m.get(), dir))
        {
            results[m[dir]] = input[i];
        }
    }
    return results;
}

mat project(const mat & input, u_int dim, u_int dir)
{
    u_int degree = 0, n_polys = input.size();
    while ((u_int) binomial_coefficient<double> (degree + dim, dim) != n_polys)
    {
        degree++;
    }
    mat results(degree + 1, vec(degree + 1, 0.));

    Multi_index_iterator m1(dim, degree), m2(dim, degree);
    u_int i,j;
    for (i = 0, m1.reset(); !m1.isFull(); i++, m1.increment())
    {
        if (isAligned(m1.get(), dir))
        {
            for (j = 0, m2.reset(); !m2.isFull(); j++, m2.increment())
            {
                if (isAligned(m2.get(), dir))
                {
                    results[m1[dir]][m2[dir]] = input[i][j];
                }
            }
        }
    }
    return results;
}

template <typename T, typename M>
T tensorize(const std::vector<M> & inputs)
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
    std::vector<T> mats(dim, eye);
    mats[dir] = matrix::convert<T,M>(input);
    return tensorize<T>(mats);
}

template std::mat tensorize(const std::vector<std::mat> & inputs);
template std::mat tensorize(const std::vector<boost::spmat> & inputs);

template boost::spmat tensorize(const std::vector<std::mat> & inputs);
template boost::spmat tensorize(const std::vector<boost::spmat> & inputs);

template std::mat tensorize(const std::mat & input, u_int dim, u_int dir);
template std::mat tensorize(const boost::spmat & input, u_int dim, u_int dir);
template boost::spmat tensorize(const std::mat & input, u_int dim, u_int dir);
template boost::spmat tensorize(const boost::spmat & input, u_int dim, u_int dir);

}
