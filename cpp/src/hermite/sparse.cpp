#include "hermite/sparse.hpp"
#include "hermite/io.hpp"

using namespace std;
namespace hermite
{

mat row_col_val(const spmat & input)
{
    mat result(3, vec(input.nnz(), 0.));
    int index = 0;
    for (cit1_t i1 = input.begin1(); i1 != input.end1(); ++i1)
    {
        for (cit2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
        {
            result[0][index] = i2.index1();
            result[1][index] = i2.index2();
            result[2][index] = *i2;
            index++;
        }
    }
    return result;
}

spmat to_spmat(const mat & input)
{
    spmat result(input.size(), input[0].size());
    for (u_int i = 0; i < input.size(); i++)
    {
        for (u_int j = 0; j < input[0].size(); j++)
        {
            if(input[i][j] != 0)
            {
                result(i, j) = input[i][j];
            }
        }
    }
    return result;
}

// spmat to_spmat(const mat & row_col_val, u_int size1, u_int size2)
// {
//     vec rows = row_col_val[0];
//     vec cols = row_col_val[1];
//     vec vals = row_col_val[2];
//     spmat result (size1, size2);
//     for (u_int i = 0; i < rows.size(); i++)
//     {
//         result(rows[i], cols[i]) = vals[i];
//     }
//     return result;
// }

spmat to_spmat(const vec & data, const ivec & indices, const ivec & indptr,
               u_int size1, u_int size2)
{
    spmat result (size1, size2);
    for (u_int i = 0; i < size1; i++)
    {
        for (u_int j = indptr[i]; j < indptr[i+1]; j++)
        {
            result(i, indices[j]) = data[j];
        }
    }
    return result;
}

mat full(const spmat & input)
{
    u_int s1 = input.size1(), s2 = input.size2();
    mat result(s1, vec(s2, 0));
    for (u_int i = 0; i < s1; i++)
    {
        for (u_int j = 0; j < s2; j++)
        {
            result[i][j] = input(i, j);
        }
    }
    return result;
}

}
