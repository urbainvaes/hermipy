#include <cmath>

#include "hermite/types.hpp"
#include "hermite/hermite.hpp"
#include "hermite/transform.hpp"
#include "hermite/iterators.hpp"
#include "hermite/varf.hpp"

#include <boost/math/special_functions/binomial.hpp>

using namespace std;

namespace hermite {

cube triple_products(u_int degree)
{
    cube products(degree + 1, mat(degree + 1, vec(2*degree + 1, 0.)));

    // Compute entries for i ≤ j ≤ k, and the rest by symmetry.
    u_int i,j,k;

    // i = 0 (orthonormality)
    for (j = 0; j <= degree; j++)
    {
        products[0][j][j] = 1.;
    }

    // i ≥ 2; j ≥ 2: h_{i+1} h_j = x REC_A(i) h_i h_j + REB_B(i) h_(i-1) h_(j)
    //                           = (REC_A(i)/REC_A(j) (h_i h_(j+1) + REC_B(j) h_i h_(j-1)) + ...
    for (i = 0; i < degree; ++i)
    {
        double ai = REC_A(i);
        double bi = REC_B(i);
        for (j = 0; j <= degree; j++)
        {
            double aj = REC_A(j);
            double bj = REC_B(j);
            for (k = 0; k <= degree; k++)
            {
                products[i+1][j][k] += j == ai/aj * products[i][j+1][k];
                products[i+1][j][k] += j == 0 ? 0. : ai/aj*bj * products[i][j-1][k];
                products[i+1][j][k] += i == 0 ? 0. : bi * products[i-1][j][k];
            }
        }
    }

    return products;
}

mat varf(
        u_int degree,
        vec const & input,
        mat const & nodes,
        mat const & weights,
        mat const & weight)
{
    u_int dim = nodes.size();

    using boost::math::binomial_coefficient;
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
    u_int n_polys_2 = (u_int) binomial_coefficient<double> (2*degree + dim, dim);

    // Products ‹h_i h_j h_k›, with (0 ≤ i,j ≤ degree) and (0 ≤ k ≤ 2*degree)
    cube products = triple_products(degree);
    cube multi_products(n_polys, mat(n_polys, vec(n_polys_2, 0.)));

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);
    Multi_index_iterator m3(dim, 2*degree);

    u_int i,j,k,l;
    for (i = 0; i <= degree; i++, m1.increment())
    {
        for (j = 0; j <= degree; j++, m2.increment())
        {
            for (k = 0; k <= degree; k++, m3.increment())
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
    vec Hf = transform(2*degree, input, nodes, weights, true, false);

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
