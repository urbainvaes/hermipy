#include <cmath>

#include <iostream>
#include <map>

#include "hermite/hermite.hpp"
#include "hermite/iterators.hpp"
#include "hermite/transform.hpp"
#include "hermite/types.hpp"
#include "hermite/varf.hpp"

#include <boost/math/special_functions/binomial.hpp>

#define MIN(i,j) (i < j ? i : j)
#define MAX(i,j) (i < j ? j : i)

using namespace std;
using boost::math::binomial_coefficient;

namespace hermite {

// Exact, does not rely on quadrature
cube triple_products_1d(int degree)
{
    cube products(degree + 1, mat(degree + 1, vec(2*degree + 1, 0.)));
    cube e_products(degree + 1, mat(2*degree + 1, vec(2*degree + 1, 0.)));

    // Compute entries for i ≤ j ≤ k, and the rest by symmetry.
    int i,j,k;

    for (j = 0; j <= 2*degree; j++)
    {
        e_products[0][j][j] = 1.;
    }

    // i ≥ 2; j ≥ 2: h_{i+1} h_j = x REC_A(i) h_i h_j - REB_B(i) h_(i-1) h_(j)
    //                           = (REC_A(i)/REC_A(j) (h_i h_(j+1) + REC_B(j) h_i h_(j-1)) - ...
    for (i = 1; i <= degree; i++)
    {
        double ai = REC_A(i-1);
        double bi = REC_B(i-1);
        for (j = i; j <= 2*degree - i; j++)
        {
            double aj = REC_A(j);
            double bj = REC_B(j);
            for (k = j - i; k <= j + i; k++)
            {
                e_products[i][j][k] += ai/aj * e_products[i-1][j+1][k];
                e_products[i][j][k] += ai/aj*bj * e_products[i-1][j-1][k];
                e_products[i][j][k] += i == 1 ? 0. : - bi * e_products[i-2][j][k];
            }
        }
    }

    for (i = 0; i <= degree; i++)
    {
        for (j = i; j <= degree; j++)
        {
            for (k = j - i; k <= j + i; k++)
            {
                products[i][j][k] = e_products[i][j][k];
                products[j][i][k] = e_products[i][j][k];
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
    u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
    u_int n_polys_2 = (u_int) binomial_coefficient<double> (2*degree + dim, dim);

    cube products = triple_products_1d(degree);

    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);
    Multi_index_iterator m3(dim, 2*degree);

    // Hermite transform of input function
    vec Hf = transform(2*degree, input, nodes, weights, true);

    u_int i,j,k,l;
    mat result(n_polys, vec(n_polys, 0.));

    for (k = 0, m3.reset(); k < n_polys_2; k++, m3.increment())
    {
        if (abs(Hf[k]) < 1e-14)
        {
            continue;
        }

        for (i = 0, m1.reset(); i < n_polys; i++, m1.increment())
        {
            for (j = 0, m2.reset(); j < n_polys; j++, m2.increment())
            {
                double increment = Hf[k];
                for (l = 0; l < dim; l++)
                {
                    increment *= products[m1[l]][m2[l]][m3[l]];
                }
                result[i][j] += increment;
            }
        }
    }

    return result;
}

string hash_print(ivec v) {
    string hash = "";
    for (u_int i = 0; i < v.size(); ++i)
    {
        hash += std::to_string(v[i]) + "-";
    }
    return hash;
}

// Only for cov = D!
mat dvarf(
        u_int dim,
        u_int degree,
        u_int direction,
        mat const & var)
{
    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);

    u_int i,j;

    map<string, u_int> lin_indices;
    for (i = 0, m1.reset(); i < var.size(); i++, m1.increment())
    {
        lin_indices.insert(pair<string, u_int>(hash_print(m1.get()), i));
    }

    mat results = mat(var.size(), vec(var.size(), 0));
    for (i = 0, m1.reset(); i < var.size(); i++, m1.increment())
    {
        for (j = 0, m2.reset(); j < var.size(); j++, m2.increment())
        {
            if (m2[direction] == 0)
            {
               continue;
            }

            ivec diff_m2 = m2.get();
            diff_m2[direction] -= 1;
            u_int id = lin_indices[hash_print(diff_m2)];
            results[i][j] = sqrt(m2[direction])*var[i][id];
            // Entry i,j correspond to < A h_j, h_i >
        }
    }
    return results;
}

}
