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

#ifndef IO_H
#define IO_H

#include "hermite/types.hpp"

#include <string>
#include <iostream>

namespace hermite {

template<typename T> 
std::ostream & printVec(std::ostream & os, const std::vector<T> & a)
{
    os << "[" << a[0];
    for (u_int i = 1; i < a.size(); i++)
    {
            os << ", " << a[i];
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> & data)
{
    printVec(os, data);
    return os;
}

template<typename T>
std::ostream &  printMat(std::ostream & os,
        const std::vector< std::vector<T> > & a,
        const std::string & begin = " ")
{
    os << "[" << a[0];
    for (u_int i = 1; i < a.size(); i++)
    {
        os << std::endl << begin << a[i];
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector< std::vector<T> > & data)
{
    printMat<T>(os, data);
    return os;
}

template<typename T> 
std::ostream &  printCube(std::ostream & os,
        const std::vector< std::vector< std::vector<T> > > & a,
        const std::string & begin = "  ")
{
    os << "[" << a[0];
    for (u_int i = 1; i < a.size(); i++)
    {
        os << std::endl << begin << a[i];
    }
    os << "]" << std::endl;
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os,
        const std::vector< std::vector< std::vector<T> >  > & data)
{
    printCube<T>(os, data);
    return os;
}

std::ostream& operator<<(std::ostream& os, const spmat & data);
std::ostream& operator<<(std::ostream& os, const boost_mat & data);

}

#endif
