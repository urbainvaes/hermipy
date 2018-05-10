// Compilation: g++ -I/usr/include/python3.6m -lpython3.6m -lboost_python3 -lboost_numpy3 test.cpp

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "hermite/sparse.hpp"

#include <iostream>

using namespace std;
using namespace hermite;
using namespace boost::numeric::ublas;
typedef compressed_matrix<double, row_major> spmat;

std::mat simple(int n)
{
    mat result(n, vec(n, 0.));
    for (u_int i = 0; i < result.size(); i++)
    {
        for (u_int j = 0; j < result[0].size(); j++)
        {
            if (i%2 == 0 && j%2 == 0)
            {
                result[i][j] = (double) i*n + (double) j;
            }
        }
    }
    return result;
}

int main()
{
    mat test = simple(6);
    spmat converted = mat_to_spmat(test);
    cout << converted << endl;
}
