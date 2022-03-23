#include <iostream>
#include <cstdlib>
#include <cmath>
#include "matrixmath2.h"


using namespace std;


int main(int argc, char *argv[]){

	int n = 3;
	double tol = 1e-10;

	double *A, *EigVal, *EigVect, *x1, *x2, *x3;
	A = new double[n*n];
	EigVal = new double[n*n];
	EigVect = new double[n*n];
	x1 = new double[3];
	x2 = new double[3];
	x3 = new double[3];

	A[0] = 1.0;
	A[1] = 1.0;
	A[2] = -2.0;

	A[3] = 0.0;
	A[4] = 0.0;
	A[5] = 2.0;

	A[6] = -1.0;
	A[7] = 0.0;
	A[8] = 1.0;

	cout << "The Matrix is:" << endl;
	print(A, n, n);

	cout << endl;
	eigenvalues(A, EigVal, n, tol, EigVect);
	
	cout<< "Its Diagonal Matrix (With Eigenvalues) is:" << endl;
	print(EigVal, n, n);
	
	cout << endl;
	cout << "Its matrix made of Eigenvectors is:" << endl;
	print(EigVect, n, n);

	for(int i=0; i<n; i++){

		x1[i] = EigVect[i+0*n];
		x2[i] = EigVect[i+1*n];
		x3[i] = EigVect[i+2*n];
	}

	cout << endl;
	cout << "When the eigenvalue is " << EigVal[0+0*n] << ", the eigenvector is:" << endl;
	print(x1, n);

	cout << endl;	
	cout << "When the eigenvalue is " << EigVal[1+1*n] << ", the eigenvector is:" << endl;
	print(x2, n);

	cout << endl;	
	cout << "When the eigenvalue is " << EigVal[2+2*n] << ", the eigenvector is:" << endl;
	print(x3, n);

	return 0;
}