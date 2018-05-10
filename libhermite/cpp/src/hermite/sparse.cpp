#include "hermite/sparse.hpp"

using namespace boost::numeric::ublas;
typedef compressed_matrix<double, row_major>::iterator1 it1_t;
typedef compressed_matrix<double, row_major>::iterator2 it2_t;
typedef compressed_matrix<double, row_major> spmat;

using namespace std;
namespace hermite
{

mat spmat_to_mat(spmat input)
{
    mat result(3, vec(input.nnz(), 0.));
    int index = 0;
    for (it1_t i1 = input.begin1(); i1 != input.end1(); ++i1)
    {
        for (it2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
        {
            result[0][index] = i2.index1();
            result[1][index] = i2.index2();
            result[2][index] = *i2;
        }
    }
    return result;
}

spmat mat_to_spmat(mat input)
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

mat full(spmat input)
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
