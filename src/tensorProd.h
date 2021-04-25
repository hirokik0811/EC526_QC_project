#ifndef TENSORPRODDENSE_H
#define TENSORPRODDENSE_H

#include <complex>
#include "sparseMat.h"


typedef std::complex<double> Complex;
typedef struct { Complex* q; int dim; } QReg;
//typedef struct { Complex* val; int* row; int* col; int nVal; int nRows; int nCols; } sparseCSR;

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


sparseCSR* tensorProdSparse(sparseCSR* a, sparseCSR* b) {
	sparseCSR* prod = new sparseCSR(a->nVal * b->nVal, a->nRows * b->nRows, a->nCols * b->nCols);
	prod->row[0] = 0;

	int p_val_ind = 0;
	for (int r_a = 0; r_a < a->nRows; r_a++) {
		for (int r_b = 0; r_b < b->nRows; r_b++) {
			for (int i = a->row[r_a]; i < a->row[r_a + 1]; i++) {
				for (int j = b->row[r_b]; j < b->row[r_b + 1]; j++) {
					prod->val[p_val_ind] = a->val[i] * b->val[j];
					prod->col[p_val_ind] = a->col[i] * b->nCols + b->col[j];
					p_val_ind++;
				}
			}
			prod->row[r_a * b->nRows + r_b + 1] = prod->row[r_a * b->nRows + r_b] + (a->row[r_a + 1] - a->row[r_a]) * (b->row[r_b + 1] - b->row[r_b]);
		}
	}

	return prod;
}

/* ATTEMPT AT A PARALLIZED VERSION -- Doesn't work due to the complex inner loops. I don't think it is possible to simplify these loop bounds without a different sparse matrix representation.

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


	#pragma acc data copy(p_val[0:p_nVal]) copy(p_row[0:p_nRows]) copy(p_col[0:p_nVal]) copyin(a_val[0:a_nVal]) copyin(a_row[0:a_nRows]) copyin(a_col[0:a_nVal]) copyin(b_val[0:b_nVal]) copyin(b_row[0:b_nRows]) copyin(b_col[0:b_nVal])
	#pragma acc parallel loop
	for (int r_a = 0; r_a < a_nRows; r_a++) {
		for (int r_b = 0; r_b < b_nRows; r_b++) {
			for (int i = a_row[r_a]; i < a_row[r_a + 1]; i++) {
				for (int j = b_row[r_b]; j < b_row[r_b + 1]; j++) {
					p_val[p_val_ind] = a_val[i] * b_val[j];
					p_col[p_val_ind] = a_col[i] * b_nCols + b_col[j];
					p_val_ind++;
				}
			}
			p_row[r_a * b_nRows + r_b + 1] = p_row[r_a * b_nRows + r_b] + (a_row[r_a + 1] - a_row[r_a]) * (b_row[r_b + 1] - b_row[r_b]);
		}
	}

	return prod;
}
*/

#endif

