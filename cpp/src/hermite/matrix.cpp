#include "hermite/matrix.hpp"

namespace hermite { namespace matrix {

    // Template specialization
    template <> mat construct(u_int size1, u_int size2)
    {
        return mat(size1, vec(size2, 0.));
    }

    template <> spmat construct(u_int size1, u_int size2)
    {
        return spmat(size1, size2);
    }

    template <> cmat construct(u_int size1, u_int size2)
    {
        auto dims = boost::extents[size1][size2];
        return cmat(dims);
    }

    template<typename T> T convert(const mat & input);
    template<typename T> T convert(const spmat & input);

    template <> mat convert (const mat & input)
    {
        return input;
    }

    template <> spmat convert (const spmat & input)
    {
        return input;
    }

    template <> mat convert(const spmat & input)
    {
        return full(input);
    }

    template <> spmat convert(const mat & input)
    {
        return to_spmat(input);
    }
}}
