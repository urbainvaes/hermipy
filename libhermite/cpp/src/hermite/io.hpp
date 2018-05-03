#ifndef IO_H
#define IO_H

#include "types.hpp"
#include <string>
#include <iostream>

namespace hermite {
using namespace std;

template<typename T> void printVec(vector<T> a)
{
    cout << "[" << a[0];
    for (u_int i = 1; i < a.size(); i++)
    {
            cout << ", " << a[i];
    }
    cout << "]";
}

template<typename T> void printMat(vector< vector<T> >  a, string begin)
{
    cout << "["; printVec<T>(a[0]);
    for (u_int i = 1; i < a.size(); i++)
    {
        cout << endl << begin;
        printVec<T>(a[i]);
    }
    cout << "]";
}

template<typename T> void printCube(vector< vector< vector<T> >  > a, string begin)
{
    cout << "[";
    printMat<T>(a[0]);
    for (u_int i = 1; i < a.size(); i++)
    {
        cout << endl << begin;
        printMat<T>(a[i]);
    }
    cout << "]" << endl;
}}


#endif
