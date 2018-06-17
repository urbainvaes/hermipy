#include "hermite/io.hpp"
#include "hermite/types.hpp"
#include "hermite/matrix.hpp"

namespace hermite
{

std::ostream& operator<<(std::ostream& os, const spmat & data)
{
    os << matrix::convert<mat>(data);
    return os;
}

std::ostream& operator<<(std::ostream& os, const boost_mat & data)
{
    os << matrix::convert<mat>(data);
    return os;
}

}
