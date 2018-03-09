#include <cmath>

#include <iostream>

#include "hermite/hermite.hpp"
#include "hermite/iterators.hpp"
#include "hermite/transform.hpp"
#include "hermite/types.hpp"
#include "hermite/varf.hpp"

#include <boost/math/special_functions/binomial.hpp>

#define MIN(i,j) (i < j ? i : j)
#define MAX(i,j) (i < j ? j : i)

using namespace std;

namespace hermite {

cube triple_products(int degree)
{
    cube products(degree + 1, mat(degree + 1, vec(2*degree + 1, 0.)));
    cube e_products(degree + 1, mat(2*degree + 1, vec(2*degree + 1, 0.)));

    // Compute entries for i ≤ j ≤ k, and the rest by symmetry.
    int i,j,k;

    // i = 0 (orthonormality)
    for (j = 0; j <= 2*degree; j++)
    {
        e_products[0][j][j] = 1.;
    }

    // i = 1, j > i
    for (j = 1; j <= 2 * degree - 1; j++)
    {
        double aj = REC_A(j);
        double bj = REC_B(j);
        e_products[1][j][j+1] = 1/aj;
        e_products[1][j][j-1] = bj/aj;
    }

    // i ≥ 2; j ≥ 2: h_{i+1} h_j = x REC_A(i) h_i h_j - REB_B(i) h_(i-1) h_(j)
    //                           = (REC_A(i)/REC_A(j) (h_i h_(j+1) + REC_B(j) h_i h_(j-1)) - ...
    for (i = 2; i <= degree; i++)
    {
        double ai = REC_A(i-1);
        double bi = REC_B(i-1);
        for (j = i; j <= 2*degree - i; j++)
        {
            double aj = REC_A(j);
            double bj = REC_B(j);
            for (k = 0; k <= i+j; k++)
            {
                e_products[i][j][k] += ai/aj * e_products[i-1][j+1][k];
                e_products[i][j][k] += ai/aj*bj * e_products[i-1][j-1][k];
                e_products[i][j][k] += - bi * e_products[i-2][j][k];
            }
        }
    }

    for (i = 0; i <= degree; i++)
    {
        for (j = 0; j <= degree; j++)
        {
            for (k = 0; k <= 2*degree; k++)
            {
                products[i][j][k] = e_products[MIN(i,j)][MAX(i,j)][k];
            }
        }
    }

    return products;
}

mat varf(
        u_int degree,
        vec const & input,
        mat const & nodes,
        mat const & weights)
{
    u_int dim = nodes.size();

    using boost::math::binomial_coefficient;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
    u_int n_polys_2 = (u_int) binomial_coefficient<double> (2*degree + dim, dim);

    // Products ‹h_i h_j h_k›, with (0 ≤ i,j ≤ degree) and (0 ≤ k ≤ 2*degree)
    cube products = triple_products((int) degree);
    cube multi_products(n_polys, mat(n_polys, vec(n_polys_2, 0.)));

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);
    Multi_index_iterator m3(dim, 2*degree);

    u_int i,j,k,l;
    for (i = 0, m1.reset(); i < n_polys; i++, m1.increment())
    {
        for (j = 0, m2.reset(); j < n_polys; j++, m2.increment())
        {
            for (k = 0, m3.reset(); k < n_polys_2; k++, m3.increment())
            {
                multi_products[i][j][k] = 1.;
                for (l = 0; l < dim; l++)
                {
                    multi_products[i][j][k] *= products[m1[l]][m2[l]][m3[l]];
                }
            }
        }
    }

    // Hermite transform of input function
    vec Hf = transform(2*degree, input, nodes, weights, true);

    mat result(n_polys, vec(n_polys, 0.));
    for (i = 0; i < n_polys; i++)
    {
        for (j = 0; j < n_polys; j++)
        {
            for (k = 0; k < n_polys_2; k++)
            {
                result[i][j] += multi_products[i][j][k] * Hf[k];
            }
        }
    }

    return result;
}

}
