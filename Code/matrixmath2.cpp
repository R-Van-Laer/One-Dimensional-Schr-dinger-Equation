#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;



void print(double *A, int nr, int nc) {
  for (int i=0; i<nr; i++) {
    for (int j=0; j<nc; j++) {
      cout << A[i + j*nr] << "\t";
    }
    cout << "\n";
  }
  cout << endl;
}

void print(double *x, int nr) {
  for (int i=0; i<nr; i++) {
      cout << x[i] << "\n";
  }
  cout << endl;
}

int gauss(double *A, int n, double *x) {

  for (int i=0; i<n; i++) {
  	// Search for maximum in this column
  	double maxl = abs(A[i + i*n]);
  	int maxr = i;
  	for (int k=i+1; k<n; k++) {
  	  if (abs(A[k + i*n]) > maxl) {
  		maxl = abs(A[k + i*n]);
  		maxr = k;
  	  }
  	}
	
  	// Swap maximum row with current row (column by column)
  	for (int k=i; k<n+1;k++) {
  	  double tmp = A[maxr + k*n];
  	  A[maxr + k*n] = A[i + k*n];
  	  A[i + k*n] = tmp;
  	}
	
  	// Make all rows below this one 0 in current column
  	for (int k=i+1; k<n; k++) {
  	  double c = -A[k + i*n]/A[i + i*n];
  	  for (int j=i; j<n+1; j++) {
    		if (i==j) {
    		  A[k + j*n] = 0;
    		} 
        else {
    		  A[k + j*n] += c * A[i + j*n];
    		}
  	  }// end of j loop
  	}// end of k loop
  }// end of i loop
  
  // Solve equation Ax=b for an upper triangular matrix A
  for (int i=n-1; i>=0; i--) {
  	x[i] = A[i + n*n]/A[i + i*n];
	 for (int k=i-1;k>=0; k--) {
	    A[k + n*n] -= A[k + i*n] * x[i];
  	}
  }

  return 0;
}

double det( double *A, int n) {
  int piv = 1; // to determine pivots
  double *Atmp;
  Atmp = new double [n*n];
  // and copy A to Atmp to not modify original matrix
  for(int i = 0 ; i < n ; i++){
 	 for(int j = 0 ; j < n ; j++){
 	 		Atmp[i + j*n] = A[i + j*n]; 
  	}
  }

  for (int i=0; i<n; i++) {
	// Search for maximum in this column
	double maxl = abs(Atmp[i + i*n]);
	int maxr = i;
	for (int k=i+1; k<n; k++) {
	  if (abs(Atmp[k + i*n]) > maxl) {
		maxl = abs(Atmp[k + i*n]);
		maxr = k;
	  }
	}
	
	// Swap maximum row with current row (column by column)
	for (int k=i; k<n+1;k++) {
	  double tmp = Atmp[maxr+ k*n];
	  Atmp[maxr+ k*n] = Atmp[i+ k*n];
	  Atmp[i+ k*n] = tmp;
	  piv *= -1;
	}
	
	// Make all rows below this one 0 in current column
	for (int k=i+1; k<n; k++) {
	  double c = -Atmp[k+ i*n]/Atmp[i+ i*n];
	  for (int j=i; j<n+1; j++) {
		if (i==j) {
		  Atmp[k+ j*n] = 0;
		} else {
		  Atmp[k+ j*n] += c * Atmp[i+ j*n];
		}
	  }
	}
  }
  
  // evaluate product of diag elements.
  double x = 1.0;
  for (int i=0; i<n; i++) {
	x*= Atmp[i+ i*n];
  }

  delete [] Atmp;


  return piv*x;
}


int matprod(double *A, double *b, double *x, int nra, int nca)  {
  // x = A * b
  for (int i=0; i<nra; i++) {
    double sum =0.0;
    for(int j = 0; j < nca; j++){
      sum += A[i + j*nra] * b[j];
    }
    x[i] = sum;
  }
  return 0;

}

int matprod(double *A, double *B, double *C, int nra, int nca, int ncb) {
  
  for (int i=0; i<nra; i++) {
	for (int k=0; k<ncb; k++) {
	  double sum =0.0;
	  for (int j=0; j<nca; j++) {
		sum += A[i + j*nra] * B[j + k*nca];
	  }
	  C[i + k*nra] = sum;
	}
  }
  return 0;

}

double inv(double *A, double *Ainv, int n) {
  int i,j,k;
  int pivotsign = 1;
  
  int *index; // to mark pivot
  double *b; 
  double *scale;
  double *Atmp;
  index = new int[n];
  scale = new double[n];
  b = new double [n*n];
  Atmp = new double [n*n];
  
 
  // Set to identity -- note new version of if/else
  // and copy A to Atmp to not modify original matrix
  for(i = 0 ; i < n ; i++){
 	 for(j = 0 ; j < n ; j++){
 	 		b[i + j*n] = (i==j) ? 1.0: 0.0; 
 	 		Atmp[i + j*n] = A[i + j*n]; 
  	}
  }

  for (int i = 0 ; i < n; i++) {
	index[i] = i; //initialize index list
	double scalemax = 0.0;
	for(j = 0 ; j < n ; j++){
	  if(scalemax > abs(Atmp[i + j*n]) )
		scalemax = scalemax;
	  else
		scalemax = abs(Atmp[i + j*n]);
	}
	scale[i] = scalemax;
  }  
  
  for (int i = 0; i < n - 1; i++) {
	// Select pivot row
	double ratiomax = 0.0;
	int jpivot = i;
	for(k = i ; k < n ; k++){
	  double ratio = abs(Atmp[index[k] + i*n]/scale[index[k]]);
	  if(ratio>ratiomax){
		jpivot = k;
		ratiomax = ratio;
	  }
	}
	// Perform pivoting using row index list
	int indexj = index[i];
	if(jpivot!=i){//PIVOT!
	  indexj = index[jpivot];
	  index[jpivot] = index[i];//Swap index jpivot and k
	  index[i] = indexj;
	  pivotsign *= -1;//flip sign of det
	}
	// Forward elimination
	for(k=i+1;k<n;k++){
	  double c = Atmp[index[k] + i*n]/Atmp[indexj + i*n];
	  for(j=i+1;j<n;j++){
		Atmp[index[k] + j*n] -= c*Atmp[indexj + j*n];
	  }
	  Atmp[index[k] + i*n] = c;
	  for(j = 0 ; j < n ; j++){
		b[index[k] + j*n] -= Atmp[index[k] + i*n]*b[indexj + j*n];
	  }
	}
  }
  double det = pivotsign;
  for(i = 0 ; i < n ; i++){
	det *= Atmp[index[i] + i*n];
  }

  // Backsubstitution
  for(k=0;k<n;k++){
	Ainv[n-1 + k*n] = b[index[n-1] + k*n]/Atmp[index[n-1] + (n-1)*n];
	for(i=n-2;i>=0;i--){
	  double sum = b[index[i] + k*n];
	  for(j=i+1;j<n;j++){
		sum -= Atmp[index[i] + j*n]*Ainv[j + k*n];
	  }
	  Ainv[i + k*n] = sum/Atmp[index[i] + i*n];
	}
  }
  delete [] Atmp;
  delete [] b;
  delete [] scale;
  delete [] index;

  return det;
}


int transpose(double *A, double *AT, int nr, int nc){
  for (int i=0; i<nr; i++) {
	for (int j=0; j<nc; j++) {
	  AT[j+ i*nc] = A[i+ j*nr];
	}
  }
  return 0;
}



int eigenvalues(double *A, double *B, int n, double tol, double *eigenvector){

	// Copy A into B;
  	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
	  		B[i + j*n] = A[i+ j*n];
		}
	}

  for(int i=0; i<n; i++){
     for(int j=0; j<n; j++){
	if(i == j){
	   eigenvector[i+j*n] = 1.0;
	}
	else{
	   eigenvector[i+j*n] = 0.0;
	}
     }
  }
	
  int count = 0, imax, jmax;
  double c, s, t, tau;
  double error = 1.0;
  while(error > tol){
	imax = 0;
	jmax = 1;
	error = abs(B[0 + 1*n]);
	for(int i = 0 ; i < n ; i++){
	  for(int j = i+1 ; j < n ; j++){
		if( abs(B[i + j*n]) > error ){
		  error = abs(B[i + j*n]);
		  imax = i;
		  jmax = j;
		}
	  }
	}	

	double diff = (B[imax + imax*n] - B[jmax  + jmax*n]);
	tau = diff/(2.0*B[imax + jmax*n]);
	t = (tau>=0) ? -tau + sqrt(1.0 + tau*tau) :-tau - sqrt(1.0 + tau*tau);

	c = 1/sqrt(1 + t*t);
	s = t*c;
	double *R, *RT, *Atmp, *Btmp;
	R = new double[n*n];
	RT = new double[n*n];
	Atmp = new double[n*n];
	Btmp = new double[n*n];

  	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
	  		R[i + j*n] = (i==j) ? 1.0 : 0.0;
		}
	}
	R[imax + imax*n] = c;
	R[jmax + jmax*n] = c;
	R[imax + jmax*n] = -s;
	R[jmax + imax*n] = s;
	transpose(R,RT,n,n);
	matprod(RT,B,Atmp,n,n,n);
	matprod(Atmp,R,B,n,n,n);
	count++ ;

	matprod(eigenvector, R, Btmp, n,n,n);

	for(int i=0; i<n; i++){
	   for(int j=0; j<n; j++){
	       eigenvector[i+j*n] = Btmp[i+j*n];
	   }
	}

	error = abs(B[0 + 1*n]);
	for(int i = 0 ; i < n ; i++){
	  for(int j = i+1 ; j < n ; j++){
		if( abs(B[i + j*n]) >= error ){
		  error = abs(B[i + j*n]);
		  imax = i;
		  jmax = j;
		}
	  }
	}
	
	error *= error;
 	delete[] Atmp;
  	delete[] RT;
  	delete[] R;
  }// end while loop
//  cout << "Took " << count << " iterations..." << endl;
  return 0;
	
}