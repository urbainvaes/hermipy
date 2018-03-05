#include <iostream>

#include "types.hpp"
#include "io.hpp"

using namespace std;

namespace hermite {

void printVec(vec a)
{
    cout << "[" << a[0];
    for (u_int i = 0; i < a.size(); i++)
    {
            cout << ", " << a[i];
    }
    cout << "]" << endl;
}

void printMat(mat a)
{
    cout << "[";
    printVec(a[0]);
    for (u_int i = 0; i < a.size(); i++)
    {
        cout << " ";
        printVec(a[i]);
    }
    cout << "]" << endl;
}

void printCube(cube a)
{
    cout << "[";
    printMat(a[0]);
    for (u_int i = 0; i < a.size(); i++)
    {
        cout << " ";
        printMat(a[i]);
    }
    cout << "]" << endl;
}

}
