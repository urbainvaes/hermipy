#ifndef LIB_H
#define LIB_H

#include <unordered_map>
#include "hermite/types.hpp"

namespace hermite
{

// Find degree associated with number of polynomials
u_int bissect_degree(u_int dim, u_int n_polys);
u_int find_dim(u_int degree, u_int n_polys);

// Hash multi-index
u_int hash_multi_ind(ivec v, int degree);
std::unordered_map<u_int, u_int> hash_table(u_int dim, u_int degree);

// Check if multi-index is aligned to another
bool isAligned(const ivec & m, u_int dir);
bool isAligned(const ivec & m, const ivec & dirs);


// Extract sub-vector
template<typename T>
std::vector<T> extract (const std::vector<T> & input, const ivec & indices)
{
    std::vector<T> result(indices.size(), 0);
    for (u_int i = 0; i < indices.size(); i++)
    {
        result[i] = input[indices[i]];
    }
    return result;
}

}

#endif
