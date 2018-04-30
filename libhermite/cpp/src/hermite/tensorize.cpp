#include <boost/math/special_functions/binomial.hpp>

#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/templates.hpp"

using namespace std;
using boost::math::binomial_coefficient;

namespace hermite {

bool isAligned(ivec m, u_int dir)
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

vec tensorize_vecs(mat input)
{
    u_int dim = input.size();
    u_int degree = input[0].size() - 1;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
    vec results(n_polys, 0.);

    Multi_index_iterator m(dim, degree);
    for (u_int i = 0; !m.isFull(); i++, m.increment())
    {
        results[i] = 1;
        for (u_int j = 0; j < dim; j++)
        {
            results[i] *= input[j][m[j]];
        }
    }
    return results;
}

// TODO: Tensorize two vectors (urbain, Sun 29 Apr 2018 09:18:08 PM BST)
vec tensorize_vec(vec input, u_int dim, u_int dir)
{
    u_int degree = input.size() - 1;
    vec cst(degree + 1, 0.); cst[0] = 1;
    mat vecs(dim, cst);
    vecs[dir] = input;
    return tensorize_vecs(vecs);
}

mat tensorize_mat(mat input, u_int dim, u_int dir)
{
    u_int degree = input.size() - 1;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
    mat results(n_polys, vec(n_polys, 0.));

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);
    u_int i,j;
    for (i = 0, m1.reset(); !m1.isFull(); i++, m1.increment())
    {
        for (j = 0, m2.reset(); !m2.isFull(); j++, m2.increment())
        {
            if (isAligned(m2.get() - m1.get(), dir))
            {
                results[i][j] = input[m1[dir]][m2[dir]];
            }
        }
    }
    return results;
}

vec project_vec(vec input, u_int dim, u_int dir)
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

mat project_mat(mat input, u_int dim, u_int dir)
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
