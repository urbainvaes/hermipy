/*
 * Copyright (C) 2018 Urbain Vaes
 *
 * This file is part of hermipy, a python/C++ library for automating the
 * Hermite Galerkin method.
 *
 * hermipy is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * hermipy is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#include "hermite/iterators.hpp"
#include "hermite/types.hpp"
#include "hermite/io.hpp"

#include <boost/math/special_functions/binomial.hpp>

using namespace hermite;
using namespace std;

int main()
{
    u_int dim = 4, degree = 100, i = 0;
    Multi_index_iterator m(dim, degree);


    for (m.reset(), i = 0; !m.isFull(); i++, m.increment())
    {
        u_int index = Multi_index_iterator::index(m.get());

        if (i != index)
        {
            return 1;
        }

        // cout << "i " << i << ", index(m_i): " << index << endl;
    }
    cout << "Test passed" << endl;
}
