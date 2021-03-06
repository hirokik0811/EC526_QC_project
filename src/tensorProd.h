#ifndef TENSORPRODDENSE_H
#define TENSORPRODDENSE_H

#include <vector>
#include <complex>
#include <Eigen/SparseCore>
//#include "sparseMat.h"

using namespace Eigen;
using namespace std;


typedef std::complex<double> Complex;
typedef Eigen::Triplet<Complex> T;
typedef Eigen::SparseMatrix<Complex,RowMajor> SpMat;
//typedef struct { Complex* q; int dim; } QReg;
//typedef struct { Complex* val; int* row; int* col; int nVal; int nRows; int nCols; } sparseCSR;

/*
QReg tensorProdDense(QReg a, QReg b) {
	// Computes the tensor product between two quantum registers.
	// A single qubit is defined as a quantum register of dimension 2.
	// Example usage:
	// QReg prod = tensorProd(&qbitA, &qbitB);
	// --> qbitA and qbitB are both of type QReg


	// Initialize output Qreg. It has a dimension of the a.dim * b.dim.
	QReg prod;
	prod.dim = a.dim * b.dim;
	prod.q = new Complex[prod.dim];

	// Store Complex array pointers and size because openACC does not handle struct objects well.
	Complex* aq = a.q;
	int aN = a.dim;

	Complex* bq = b.q;
	int bN = b.dim;

	Complex* prodq = prod.q;
	int prodN = prod.dim;

	// Perform tensor product to populate prod
#pragma acc data copy(prodq[0:prodN]) copyin(aq[0:aN]) copyin(bq[0:bN])
#pragma acc parallel loop
	for (int i = 0; i < aN; i++) {
		for (int j = 0; j < bN; j++) {
			prodq[i * bN + j] = aq[i] * bq[j];
		}
	}

	return prod;
}
*/


SpMat tensorProdSparse(SpMat A, SpMat B) {
	std::vector<T> tripletList;
	tripletList.reserve(A.nonZeros() * B.nonZeros());

	for (int r_A = 0; r_A < A.outerSize(); ++r_A) {
		for (int r_B = 0; r_B < B.outerSize(); ++r_B) {
			for (SpMat::InnerIterator it_A(A, r_A); it_A; ++it_A) {
				for (SpMat::InnerIterator it_B(B, r_B); it_B; ++it_B) {
					tripletList.push_back(T(it_A.row() * B.outerSize() + it_B.row(), it_A.col() * B.innerSize() + it_B.col(), it_A.value() * it_B.value()));
				}
			}
		}
	}

	SpMat prod(A.rows() * B.rows(), A.cols() * B.cols());
	prod.setFromTriplets(tripletList.begin(), tripletList.end());
	return prod;
}


/*
sparseCSR* tensorProdSparse(sparseCSR* a, sparseCSR* b) {

	int a_nVal = a->nVal;
	int a_nRows = a->nRows;
	int a_nCols = a->nCols;
	Complex* a_val = a->val;
	int* a_row = a->row;
	int* a_col = a->col;

	int b_nVal = b->nVal;
	int b_nRows = b->nRows;
	int b_nCols = b->nCols;
	Complex* b_val = b->val;
	int* b_row = b->row;
	int* b_col = b->col;


	int p_nVal = a_nVal * b_nVal;
	int p_nRows = a_nRows * b_nRows;
	int p_nCols = a_nCols * b_nCols;
	sparseCSR* prod = new sparseCSR(p_nVal, p_nRows, p_nCols);
	Complex* p_val = prod->val;
	int* p_row = prod->row;
	int* p_col = prod->col;
	p_row[0] = 0;



	int p_val_ind = 0;
#pragma acc data copyout(p_val[0:p_nVal], p_row[0:p_nRows+1], p_col[0:p_nVal]) copyin(a_val[0:a_nVal], a_row[0:a_nRows+1], a_col[0:a_nVal], b_val[0:b_nVal], b_row[0:b_nRows+1], b_col[0:b_nVal], p_val_ind, a_nRows, b_nRows, b_nCols)
#pragma acc region
	{

#pragma acc loop independent seq
		for (int r_a = 0; r_a < a_nRows; r_a++) {
#pragma acc loop independent seq
			for (int r_b = 0; r_b < b_nRows; r_b++) {
#pragma acc loop independent seq
				for (int i = a_row[r_a]; i < a_row[r_a + 1]; i++) {
#pragma acc loop independent seq
					for (int j = b_row[r_b]; j < b_row[r_b + 1]; j++) {
						p_val[p_val_ind] = a_val[i] * b_val[j];
						p_col[p_val_ind] = a_col[i] * b_nCols + b_col[j];
						p_val_ind++;
					}
				}
				p_row[r_a * b_nRows + r_b + 1] = p_row[r_a * b_nRows + r_b] + (a_row[r_a + 1] - a_row[r_a]) * (b_row[r_b + 1] - b_row[r_b]);
			}
		}
	}

	return prod;
}
*/

#endif

