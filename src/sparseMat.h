#ifndef SPARSEMAT_H
#define SPARSEMAT_H

#include <complex>

typedef std::complex<double> Complex;
class sparseCSR { 
public: 
  Complex* val; int* row; int* col; int nVal; int nRows; int nCols; 
  sparseCSR() {} // empty object. user must allocate explicitly. 
  sparseCSR(int nvals, int dim) {  // square matrix
    nVal = nvals; 
    nRows = dim; 
    nCols = dim; 
    val = new Complex[nVal]; 
    row = new int[dim+1];
    col = new int[nVal];
  }
  sparseCSR(int nvals, int nrows, int ncols) {  // non-square matrix
    nVal = nvals; 
    nRows = nrows; 
    nCols = nCols; 
    val = new Complex[nVal]; 
    row = new int[nRows+1];
    col = new int[nVal];
  }
  ~sparseCSR() {
    delete[] val; delete[] row; delete[] col; 
  }
  sparseCSR copy() {
    sparseCSR obj(nVal, nRows, nCols); 
    std::copy(val, val + nVal, obj.val); 
    std::copy(row, row + nRows+1, obj.row); 
    std::copy(col, col + nVal, obj.col); 
    return obj; 
  }
};

#endif
