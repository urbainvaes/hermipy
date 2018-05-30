#include "hermite/iterators.hpp"
#include "hermite/types.hpp"
#include "hermite/io.hpp"

#include <boost/math/special_functions/binomial.hpp>

using namespace hermite;
using namespace std;

int main()
{
    u_int dim = 4, degree = 10;
    Multi_index_iterator m(dim, degree);

    for (u_int i = 0; !m.isFull(); i++, m.increment())
    {
        u_int index = Multi_index_iterator::index(m.get());
        if (i != index)
        {
            return 1;
        }

        cout << "i " << i << ", index(m_i): " << index << endl;
    }
    cout << "Test passed" << endl;
}
