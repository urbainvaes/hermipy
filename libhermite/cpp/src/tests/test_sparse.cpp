#include <iostream>
#include <fstream>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

int main()
{
    using std::cerr;
    using std::cout;
    using std::endl;
    using namespace boost::numeric::ublas;
    typedef compressed_matrix<float, row_major> cMatrix;

    const size_t size = 5;
    const size_t rowInd[5] = { 0, 0, 1, 2, 4 };
    const size_t colInd[5] = { 0, 2, 0, 4, 4 };

    cMatrix sparseMat(size,size);
    for (size_t i=0; i<size; ++i)
        sparseMat(rowInd[i], colInd[i]) = 0.0;

    cout << sparseMat(0, 1) << endl;
    cout << sparseMat << endl;

    // Try writing to file
    std::ofstream ofsDat("temp.dat", std::ios::out | std::ios::binary);
    for(cMatrix::const_iterator1 rowIter = sparseMat.begin1(); rowIter != sparseMat.end1(); ++rowIter)  {
        for(cMatrix::const_iterator2 colIter = rowIter.begin(); colIter != rowIter.end(); ++colIter)    {
            ofsDat << " " << colIter.index1() << " " << colIter.index2() << " " << *colIter;
        }       // end for colIter
    }       // end for rowIter
    ofsDat.close();

    cout << "Writing ended, starting to read" << endl;

    // Try reading the file
    cMatrix sparseMat_2(size, size);
    std::ifstream ifsDat("temp.dat", std::ios::in | std::ios::binary);
    size_t rowTemp, colTemp;
    float valTemp;
    while(!ifsDat.eof())    {
        ifsDat >> rowTemp >> colTemp >> valTemp;
        cout << "row " << rowTemp << " column " << colTemp << " value " << valTemp << endl;
        sparseMat_2.insert_element(rowTemp, colTemp, valTemp);
    }

    cout << sparseMat_2 << endl;

    cout << "Test passed" << endl;
    return 0;
}
