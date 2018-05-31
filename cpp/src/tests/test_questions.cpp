#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <iostream>

typedef std::vector<double> vec;
typedef std::vector<vec> mat;
typedef boost::numeric::ublas::compressed_matrix<double, boost::numeric::ublas::row_major> spmat;

namespace matrix
{
    inline void set(spmat & input, u_int i, u_int j, double val)
    {
        input(i, j) = val;
    }

    inline void set(mat & input, u_int i, u_int j, double val)
    {
        input[i][j] = val;
    }

    inline u_int size1(const mat & input)
    {
        return input.size();
    }

    inline u_int size2(const mat & input)
    {
        return input[0].size();
    }

    inline u_int size1(const spmat & input)
    {
        return input.size1();
    }

    inline u_int size2(const spmat & input)
    {
        return input.size2();
    }

    inline double get(const spmat & input, u_int i, u_int j)
    {
        return input(i, j);
    }

    inline double get(const mat & input, u_int i, u_int j)
    {
        return input[i][j];
    }
}

int main()
{
    std::cout << "Hello, world!" << std::endl;
}
