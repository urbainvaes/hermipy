#ifndef IO_H
#define IO_H

#include "types.hpp"
#include <string>
#include <iostream>

namespace hermite {

template<typename T> std::ostream & printVec(std::ostream & os,
        std::vector<T> a)
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
std::ostream& operator<<(std::ostream& os,
        const std::vector<T> & data)
{
    printVec(os, data);
    return os;
}

template<typename T> std::ostream &  printMat(std::ostream & os,
        std::vector< std::vector<T> >  a,
        std::string begin = " ")
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
std::ostream& operator<<(std::ostream& os,
        const std::vector< std::vector<T> > & data)
{
    printMat<T>(os, data);
    return os;
}

template<typename T> std::ostream &  printCube(std::ostream & os,
        std::vector< std::vector< std::vector<T> >  > a,
        std::string begin = "  ")
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

}


#endif
