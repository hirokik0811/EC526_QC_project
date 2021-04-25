#include <iostream>
#include "tensorProd.h"

using namespace std;

int main()
{
	sparseCSR A(1, 2, 1);
	sparseCSR B(1, 2, 1);
	sparseCSR* C;

	// A = |0>
	A.val[0] = 1;
	A.row[0] = 0; A.row[1] = 1; A.row[2] = 1;
	A.col[0] = 0;

	// B = |1>
	B.val[0] = 1;
	B.row[0] = 0; B.row[1] = 0; B.row[2] = 1;
	B.col[0] = 0;

	#pragma acc init

	// C = |01>
	C = tensorProdSparse(&A, &B);

	
	C->print();

	delete C;

	return 0;
}
