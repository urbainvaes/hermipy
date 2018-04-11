#include <iostream>

#include "types.hpp"
#include "io.hpp"

using namespace std;

namespace hermite {

void printVec(vec a)
{
    cout << "[" << a[0];
    for (u_int i = 1; i < a.size(); i++)
    {
            cout << ", " << a[i];
    }
    cout << "]";
}

void printMat(mat a, string begin)
{
    cout << "["; printVec(a[0]);
    for (u_int i = 1; i < a.size(); i++)
    {
        cout << endl << begin;
        printVec(a[i]);
    }
    cout << "]";
}

void printCube(cube a, string begin)
{
    cout << "[";
    printMat(a[0]);
    for (u_int i = 1; i < a.size(); i++)
    {
        cout << endl << begin;
        printMat(a[i]);
    }
    cout << "]" << endl;
}}
