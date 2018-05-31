#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "hermite/types.hpp"
#include "hermite/sparse.hpp"

namespace hermite { 
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

    template<typename T, typename S> T convert(const S & input);
    template<typename T> T construct(u_int size1, u_int size2);
    template<typename T> T eye(u_int size)
    {
        T result = construct<T> (size, size);
        for (u_int i = 0; i < size; i++)
        {
            set(result, i, i, 1.);
        }
        return result;
    }
}}

#endif
