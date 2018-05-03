#include <cmath>

#include <iostream>
#include <unordered_map>
#include <ctime>

// #include "sparse-matrix/SparseMatrix.h"

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

// SparseMatrix<double> sparseVarf(
//         u_int degree,
//         vec const & input,
//         mat const & nodes,
//         mat const & weights)
// {
    // u_int dim = nodes.size();
    // u_int n_polys = (u_int) binomial_coefficient<double> (degree + dim, dim);
    // u_int n_polys_2 = (u_int) binomial_coefficient<double> (2*degree + dim, dim);

    // cube products = triple_products_1d(degree);

    // Multi_index_iterator m1(dim, degree);
    // Multi_index_iterator m2(dim, degree);
    // Multi_index_iterator m3(dim, 2*degree);

    // // Hermite transform of input function
    // vec Hf = transform(2*degree, input, nodes, weights, true);

    // u_int i,j,k,l;
    // SparseMatrix<double> result(n_polys, n_polys);

    // for (k = 0, m3.reset(); k < n_polys_2; k++, m3.increment())
    // {
    //     if (abs(Hf[k]) < 1e-14)
    //     {
    //         continue;
    //     }

    //     for (i = 0, m1.reset(); i < n_polys; i++, m1.increment())
    //     {
    //         for (j = 0, m2.reset(); j < n_polys; j++, m2.increment())
    //         {
    //             double increment = Hf[k];
    //             for (l = 0; l < dim; l++)
    //             {
    //                 increment *= products[m1[l]][m2[l]][m3[l]];
    //             }
    //             result.set(result.get(i+1, j+1) + increment, i+1, j+1);
    //         }
    //     }
    // }

    // return result;
// }

c_mat contig_mat(int rows, int cols)
{
    auto dims = boost::extents[rows][cols];
    return boost::multi_array<double, 2>(dims);
}

c_mat varf(
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
    c_mat result = contig_mat(n_polys, n_polys);

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

u_int hash(ivec v, int degree) {
    u_int base = degree + 1;
    u_int result = 0;
    u_int unit = 1;
    for(u_int i = 0; i < v.size(); i++)
    {
        result += v[i]*unit;
        unit *= base;
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
mat varfd(
        u_int dim,
        u_int degree,
        u_int direction,
        mat const & var)
{
    Multi_index_iterator m1(dim, degree);
    Multi_index_iterator m2(dim, degree);

    u_int i,j;

    unordered_map<u_int, u_int> lin_indices;
    for (i = 0, m1.reset(); i < var.size(); i++, m1.increment())
    {
        lin_indices.insert(pair<u_int, u_int>(hash(m1.get(), degree), i));
    }

    mat results = mat(var.size(), vec(var.size(), 0));
    for (j = 0, m2.reset(); j < var.size(); j++, m2.increment())
    {
        if (m2[direction] == 0)
        {
           continue;
        }

        ivec diff_m2 = m2.get();
        diff_m2[direction] -= 1;
        u_int id = lin_indices[hash(diff_m2, degree)];

        for (i = 0, m1.reset(); i < var.size(); i++, m1.increment())
        {
            results[i][j] = var[i][id]*sqrt(m2[direction]);
            // Entry i,j correspond to < A h_j, h_i >
        }
    }

    return results;
}

}
