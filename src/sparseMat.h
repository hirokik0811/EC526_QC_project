#ifndef SPARSEMAT_H
#define SPARSEMAT_H

#include <iostream>
#include <complex>
#include <algorithm>

using namespace std;

typedef std::complex<double> Complex;
class sparseCSR {
public:
    Complex* val; int* row; int* col; int nVal; int nRows; int nCols;
    sparseCSR() { } // empty object. user must allocate explicitly. 
    sparseCSR(int nvals, int dim) {  // square matrix
        nVal = nvals;
        nRows = dim;
        nCols = dim;
        val = new Complex[nVal];
        row = new int[dim + 1];
        col = new int[nVal];
    }
    sparseCSR(int nvals, int nrows, int ncols) {  // non-square matrix
        nVal = nvals;
        nRows = nrows;
        nCols = nCols;
        val = new Complex[nVal];
        row = new int[nRows + 1];
        col = new int[nVal];
    }

    ~sparseCSR() {
        delete[] val; delete[] row; delete[] col;
    }

    /*
    sparseCSR copy() {
        sparseCSR obj(nVal, nRows, nCols);
        copy(val, val + nVal, obj.val);
        copy(row, row + nRows + 1, obj.row);
        copy(col, col + nVal, obj.col);
        return obj;
    }
    */

    sparseCSR(const sparseCSR& obj) {
        // Copy constructor
        nVal = obj.nVal;
        nRows = obj.nRows;
        nCols = obj.nCols;
        val = new Complex[obj.nVal];
        row = new int[obj.nRows + 1];
        col = new int[obj.nVal];

        copy(obj.val, obj.val + obj.nVal, val);
        copy(obj.row, obj.row + obj.nRows + 1, row);
        copy(obj.col, obj.col + obj.nVal, col);
    }

    void print() {
        // Prints the val, row, and col arrays
        cout << "val = ";
        for (int i = 0; i < nVal; i++) cout << val[i] << ", ";
        cout << "\n";

        cout << "row = ";
        for (int i = 0; i < nRows + 1; i++) cout << row[i] << ", ";
        cout << "\n";

        cout << "col = ";
        for (int i = 0; i < nVal; i++) cout << col[i] << ", ";
        cout << "\n";
    }
};

#endif
