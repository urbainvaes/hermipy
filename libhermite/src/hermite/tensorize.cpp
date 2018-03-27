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

vec tensorize_vec(vec input, u_int dim, u_int dir)
{
    u_int degree = input.size() - 1;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
    vec results(n_polys, 0.);

    Multi_index_iterator m(dim, degree);
    for (u_int i = 0; !m.isFull(); i++, m.increment())
    {
        if (isAligned(m.get(), dir))
        {
            results[i] = input[m[dir]];
        }
    }
    return results;
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

}
