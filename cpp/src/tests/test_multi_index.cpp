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

using namespace hermite;
using namespace std;

int main()
{
    u_int dim = 2, degree = 9, i = 0;

    Rectangle_iterator m_rect(dim, degree);
    for (m_rect.reset(), i = 0; !m_rect.isFull(); i++, m_rect.increment())
    {
        u_int index = m_rect.index(m_rect.get());
        cout << "i " << i << ", index(m_rect): " << index << ", m_rect: " << m_rect.get() << endl;
        if (i != index)
            return 1;
    }

    Cross_iterator mh(dim, degree);
    for (mh.reset(), i = 0; !mh.isFull(); i++, mh.increment())
    {
        u_int index = mh.index(mh.get());
        cout << "i " << i << ", index(mh): " << index << ", mh: " << mh.get() << endl;
        if (i != index)
            return 1;
    }

    Cross_iterator_nc mh_nc(dim, degree);
    for (mh_nc.reset(), i = 0; !mh_nc.isFull(); i++, mh_nc.increment())
    {
        u_int index = mh_nc.index(mh_nc.get());
        cout << "i " << i << ", index(mh_nc): " << index << ", mh_nc: " << mh_nc.get() << endl;
        if (i != index)
            return 1;
    }

    Triangle_iterator m(dim, degree);
    if (m.s_size(dim, degree) != m.size())
    {
        return 1;
    }

    for (m.reset(), i = 0; !m.isFull(); i++, m.increment())
    {
        u_int index = m.s_index(m.get());
        cout << "i " << i << ", index(mt): " << index << ", mt: " << m.get() << endl;
        if (i != index)
            return 1;
    }

    Cube_iterator mc(dim, degree);
    if (mc.s_size(dim, degree) != mc.size())
    {
        return 1;
    }

    for (mc.reset(), i = 0; !mc.isFull(); i++, mc.increment())
    {
        u_int index = mc.s_index(mc.get());
        cout << "i " << i << ", index(mc): " << index << ", mc: " << mc.get() << endl;
        if (i != index)
            return 1;
    }
    cout << "Test passed" << endl;
}
