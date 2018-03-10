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

mat varf_axes(
        int degree, // ! 0 - (unsigned int) 1 = +inf
        vec const & input,
        mat const & nodes,
        mat const & weights)
{
    mat products(degree + 1, vec(1*degree + 1, 0.));
    mat e_products(degree + 1, vec(2*degree + 1, 0.));

    // Compute entries for i ≤ j, and the rest by symmetry.
    int i,j;

    // i = 0: Fast Hermite transformunction
    vec Hf = transform(2*degree, input, nodes, weights, true);

    // i ≥ 2; j ≥ 2: h_{i+1} h_j = x REC_A(i) h_i h_j - REB_B(i) h_(i-1) h_(j)
    //                           = (REC_A(i)/REC_A(j) (h_i h_(j+1) + REC_B(j) h_i h_(j-1)) - ...
    for (i = 0; i <= degree; i++)
    {
        double ai = REC_A(i-1);
        double bi = REC_B(i-1);
        for (j = i; j <= 2*degree - i; j++)
        {
            if(i == 0)
            {
                if (abs(Hf[j]) > 1e-14)
                {
                    e_products[i][j] = Hf[j];
                }
                else
                {
                    e_products[i][j] = 0.;
                }
            }
            else
            {
                double aj = REC_A(j);
                double bj = REC_B(j);
                e_products[i][j] += ai/aj * e_products[i-1][j+1];
                e_products[i][j] += ai/aj*bj * e_products[i-1][j-1];
                e_products[i][j] += i == 1 ? 0. : - bi * e_products[i-2][j];
            }
        }
    }

    for (i = 0; i <= degree; i++)
    {
        for (j = 0; j <= degree; j++)
        {
            products[i][j] = e_products[MIN(i,j)][MAX(i,j)];
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

    mat products = varf_axes(degree, input, nodes, weights);
    mat result(n_polys, vec(n_polys, 0.));

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);

    u_int i,j,k;
    for (i = 0, m1.reset(); i < n_polys; i++, m1.increment())
    {
        for (j = 0, m2.reset(); j < n_polys; j++, m2.increment())
        {
            result[i][j] = 1.;
            for (k = 0; k < dim; k++)
            {
                result[i][j] *= products[m1[k]][m2[k]];
            }
        }
    }

    return result;
}

}
