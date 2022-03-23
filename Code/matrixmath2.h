#ifndef MATRICES_H
#define MATRICES_H

#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;


void print(double *A, int nr, int nc);
void print(double *x, int nr);

int gauss(double *A, int n, double *x);
double det( double *A, int n);

int matprod(double *A, double *b, double *x, 
	int nra, int nca) ;

int matprod(double *A, double *B, double *C, 
	int nra, int nca, int ncb) ;
double inv(double *A, double *Ainv, int n);
int transpose(double *A, double *AT, int nr, int nc);
int eigenvalues(double *A, double *B, int n, double tol, double *eigenvector);


#endif
