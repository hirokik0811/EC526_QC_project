#include <iostream>
#include <vector>
#include "tensorProd.h"

using namespace std;

int main()
{
	SpMat A(2, 1);
	SpMat B(2, 1);

	SpMat C;

	std::vector<T> tripletList;
	//tripletList.reserve(1);

	// A = |0>
	tripletList.push_back(T(0,0,1));
	A.setFromTriplets(tripletList.begin(), tripletList.end());
	A.makeCompressed();
	tripletList.clear();

	// B = |1>
	tripletList.push_back(T(1,0,1));
	B.setFromTriplets(tripletList.begin(), tripletList.end());

	//#pragma acc init

	// C = |01>
	C = tensorProdSparse(A, B);

	cout << C << endl;

	return 0;
}
