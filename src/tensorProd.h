#ifndef TENSORPRODDENSE_H
#define TENSORPRODDENSE_H

#include <iostream>
#include <complex>

typedef std::complex<double> Complex;
typedef struct { Complex* q; int dim; } QReg;
typedef struct { Complex* val; int* row; int* col; int nVal; int nRows; int nCols; } sparseCSR;

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


sparseCSR tensorProdSparse(sparseCSR a, sparseCSR b) {
	sparseCSR prod;

	prod.nVal = a.nVal * b.nVal;
	prod.nRows = a.nRows * b.nRows;
	prod.nCols = a.nCols * b.nCols;
	prod.val = new Complex[prod.nVal];
	prod.row = new int[prod.nRows];
	prod.col = new int[prod.nVal];

	prod.row[0] = 0;

	int row_length_a,
		row_length_b,
		a_val_row_start = 0,
		b_val_row_start = 0,
		p_val_ind = 0;

	for (int r_a = 0; r_a < a.nRows; r_a++) {
		row_length_a = a.row[r_a + 1] - a.row[r_a];
		b_val_row_start = 0;

		for (int r_b = 0; r_b < b.nRows; r_b++) {
			row_length_b = b.row[r_b + 1] - b.row[r_b];

			for (int i = 0; i < row_length_a; i++) {
				for (int j = 0; j < row_length_b; j++) {
					prod.val[p_val_ind] = a.val[a_val_row_start + i] * b.val[b_val_row_start + j];
					prod.col[p_val_ind] = a.col[a_val_row_start + i] * b.nCols + b.col[b_val_row_start + j];
					p_val_ind++;
				}
			}
			prod.row[r_a * b.nRows + r_b + 1] = prod.row[r_a * b.nRows + r_b] + row_length_a * row_length_b;
			b_val_row_start += row_length_b;
		}
		a_val_row_start += row_length_a;
	}

	return prod;
}


#endif
