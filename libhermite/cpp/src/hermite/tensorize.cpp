#include <boost/math/special_functions/binomial.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <iostream>

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

boost::spmat tensorize(const std::vector<boost::spmat> & inputs)
{
    u_int dim = inputs.size();

    if (dim == 1)
    {
        return inputs[0];
    }

    u_int degree = inputs[0].size1() - 1;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);

    boost::spmat product(n_polys, n_polys);

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
                 result *= inputs[k](m1[k], m2[k]);
            }
            if (result != 0)
            {
                product(i,j) = result;
            }
        }
    }
    return product;
}

mat tensorize(const cube & inputs)
{
    u_int dim = inputs.size();

    if (dim == 1)
    {
        return inputs[0];
    }

    u_int degree = inputs[0].size() - 1;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
    mat results(n_polys, vec(n_polys, 0.));

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);
    u_int i,j,k;
    for (i = 0, m1.reset(); !m1.isFull(); i++, m1.increment())
    {
        for (j = 0, m2.reset(); !m2.isFull(); j++, m2.increment())
        {
            results[i][j] = 1.;
            for (k = 0; k < dim; k++)
            {
                results[i][j] *= inputs[k][m1[k]][m2[k]];
            }
        }
    }
    return results;
}

mat tensorize(const mat & input, u_int dim, u_int dir)
{
    u_int degree = input.size() - 1;
    vec zeros(degree + 1, 0.);
    mat eye(degree + 1, zeros);
    for (u_int i = 0; i < degree + 1; ++i)
    {
        eye[i][i] = 1.;
    }
    cube mats(dim, eye);
    mats[dir] = input;
    return tensorize(mats);
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

}
