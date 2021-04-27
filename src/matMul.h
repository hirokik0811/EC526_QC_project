#ifndef MATRIXMUL_H
#define MATRIXMUL_H

#include "sparseMat.h"

SpMat matrixMul(SpMat A, SpMat B)
{
    int m = A.rows();
    int n = B.cols();
    SpMat C(m, n);
    C = A * B;
    return C;
}

#endif