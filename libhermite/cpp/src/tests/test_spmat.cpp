// Compilation: g++ -I/usr/include/python3.6m -lpython3.6m -lboost_python3 -lboost_numpy3 test.cpp

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "hermite/sparse.hpp"
#include "hermite/iterators.hpp"
#include "hermite/tensorize.hpp"
#include "hermite/templates.hpp"
#include "hermite/io.hpp"

#include <iostream>

using namespace std;
using namespace hermite;
namespace bn = boost::numeric::ublas;
typedef bn::compressed_matrix<double, bn::row_major> spmat;

std::mat simple(int n)
{
    mat result(n, vec(n, 0.));
    for (u_int i = 0; i < result.size(); i++)
    {
        for (u_int j = 0; j < result[0].size(); j++)
        {
            if (true || i%2 == 0 && j%2 == 0)
            {
                result[i][j] = (double) i*n + (double) j + 1;
            }
        }
    }
    return result;
}

int main()
{
    int degree = 3;
    mat initial = simple(degree + 1);
    cout << initial << endl;
    spmat sp_initial = mat_to_spmat(initial);
    std::cube inputs(2, initial);
    std::vector<spmat> sp_inputs(2, sp_initial);
    spmat sp_tensorized = sp_tensorize_mats(sp_inputs);
    mat tensorized = tensorize_mats(inputs);
    mat difference = hermite::full(sp_tensorized) - tensorized;
    imat m = list_multi_indices(2, degree);
    for (u_int i = 0; i < tensorized.size(); i++)
    {
        for (u_int j = 0; j < tensorized.size(); j++)
        {
            cout << m[i] << ", " << m[j] << ", "
                << initial[m[i][0]][m[j][0]] <<  " * "
                << initial[m[i][1]][m[j][1]] << " = "
                <<  tensorized[i][j] << "?" << endl;
        }
    }
}
