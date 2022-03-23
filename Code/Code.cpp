#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "matrixmath2.h"
#include "Potential.h"


using namespace std;


double create_mat(double L, double R, int n, double *Ham, double *x_var,
		   double (*Pot)(double x, double *par),
		   double *par){

	double tau, h2over2m;
	tau = abs(R - L) / ( (double) n);
	h2over2m = 1.0;

	for(int i=0; i<n; i++){

		x_var[i] = L + i*tau;
	}
cout << "x_var done, x_var[n] = " << x_var[n] << endl;
	for(int i=0; i<n; i++){

		for(int j=0; j<n; j++){

			if( i==j ){

				Ham[i + j*n] = (-2.0 * -h2over2m) + Pot(x_var[i], par);
			}
			else if( i==(j+1) || j==(i+1) ){

				Ham[i + j*n] = (1.0 * -h2over2m);
			}
			else{
				
				Ham[i + j*n] = 0.0;
			}
		}
	}
cout << "Ham done, Ham[1,1] = " << Ham[0] << endl;	
	return 0;
}


int main(int argc, char *argv[]){


	int Type;
	int n;
	double L, R;

	double *par;
	par = new double[3];
	for(int i=0; i<3; i++){

		par[i] = 0.0;
	}

	string Fileps, Filedat, File;

	double *x_var;
	double *Ham;
	double tol = 1e-4;
	double *eigenvalue, *eigenvector;


	if(argc < 3){

		cout << "Usage: " << argv[0] << " Type n " << endl;
		cout << "n is the size of the data set." << endl;
		cout << "This code will graph the wavefunctions with a dataset of nMAX = 100" << endl;
		cout << "Type is a integer input corresponding to one potential function." << endl;
		cout << "\t" << "1 ----- Quantum Harmonic Oscillator" << endl;
		cout << "\t" << "2 ----- Finite Potential Well" << endl;
		cout << "\t" << "3 ----- Free Particle" << endl;
		cout << "\t" << "4 ----- Constant Potential" << endl;
		cout << "\t" << "5 ----- Tunneling" << endl;
		cout << "\t" << "6 ----- Potential Step" << endl;
		return -1;
	}

	else{

		Type = atoi(argv[1]);
		n = atoi(argv[2]);
	}


	if(Type == 1){

		L = -1.0;
		R = 1.0;

		par[0] += 1.0;
		par[1] += 1.0;

		x_var = new double[n];
		Ham = new double[n*n];

		eigenvalue = new double[n*n];
		eigenvector = new double[n*n];

		create_mat(L, R, n, Ham, x_var, QHO, par);

		eigenvalues(Ham, eigenvalue, n, tol, eigenvector);

		cout << "Eigenvalues done" << endl;
	}

	if(Type == 2){

		L = 0.0;
		R = 3.0;

		par[0] += 1.0;
		par[1] += 2.0;
		par[2] += 50.0;

		x_var = new double[n];
		Ham = new double[n*n];

		eigenvalue = new double[n*n];
		eigenvector = new double[n*n];

		create_mat(L, R, n, Ham, x_var, Well, par);

		eigenvalues(Ham, eigenvalue, n, tol, eigenvector);

		cout << "Eigenvalues done" << endl;
	}

	if(Type == 3){

		L = 0.0;
		R = 3.0;

		x_var = new double[n];
		Ham = new double[n*n];

		eigenvalue = new double[n*n];
		eigenvector = new double[n*n];

		create_mat(L, R, n, Ham, x_var, Free, par);

		eigenvalues(Ham, eigenvalue, n, tol, eigenvector);
		
		cout << "Eigenvalues done" << endl;
		//print(eigenvalue, n, n);
		//print(eigenvector, n, n);
	}

	if(Type == 4){

		L = 0.0;
		R = 3.0;
		
		par[0] += 50.0;

		eigenvalue = new double[n*n];
		eigenvector = new double[n*n];

		x_var = new double[n];
		Ham = new double[n*n];

		create_mat(L, R, n, Ham, x_var, Const, par);

		eigenvalues(Ham, eigenvalue, n, tol, eigenvector);

		cout << "Eigenvalues done" << endl;
	}

	if(Type == 5){

		L = 0.0;
		R = 3.0;

		par[0] += 1.0;
		par[1] += 2.0;
		par[2] += 50.0;

		x_var = new double[n];
		Ham = new double[n*n];

		eigenvalue = new double[n*n];
		eigenvector = new double[n*n];

		create_mat(L, R, n, Ham, x_var, Tunn, par);

		eigenvalues(Ham, eigenvalue, n, tol, eigenvector);

		cout << "Eigenvalues done" << endl;
	}

	if(Type == 6){

		L = 0.0;
		R = 3.0;

		par[0] = 1.5;
		par[1] = 50.0;

		x_var = new double[n];
		Ham = new double[n*n];

		eigenvalue = new double[n*n];
		eigenvector = new double[n*n];

		create_mat(L, R, n, Ham, x_var, Step, par);

		eigenvalues(Ham, eigenvalue, n, tol, eigenvector);

		cout << "Eigenvalues done" << endl;
	}


	double sum, norm;

	for(int j=0; j<n; j++){

		sum = 0.0;
		norm = 0.0;
		
		for(int i=0; i<n; i++){

			sum += eigenvector[i + j*n] * eigenvector[i + j*n];
		}

		norm += sqrt(sum);
		
		for(int i=0; i<n; i++){

			eigenvector[i + j*n] = eigenvector[i + j*n] / norm;
		}

		norm = 0.0;
		sum = 0.0;

		cout << "Normalization done" << endl;
	}


	Filedat = File + ".dat";
	Fileps = File + ".ps";
	
	ofstream fout("DatFile.dat");

	for(int i=0; i<n; i++){
	
		for(int j=0; j<n; j++){

			fout << x_var[i] << "\t" << eigenvector[i + j*n];
		}
		
		fout << endl;
	}

	cout << "DATA File created" << endl;
	
	fout.close();


	FILE *gnuplot = popen("gnuplot", "w");
	fprintf(gnuplot, "set out 'DatFile.ps'\n");
	fprintf(gnuplot, "set term post land\n");
	fprintf(gnuplot, "set size square\n");
	fprintf(gnuplot, "set title 'Psi vs x'\n");
	fprintf(gnuplot, "set xlabel 'x'\n");
	fprintf(gnuplot, "set ylabel 'Psi(x)'\n");

	fprintf(gnuplot, "plot 'DatFile.dat' u 1:2 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:3 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:4 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:5 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:6 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:7 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:8 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:9 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:10 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:11 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:12 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:13 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:14 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:15 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:16 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:17 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:18 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:19 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:20 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:21 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:22 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:23 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:24 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:25 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:26 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:27 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:28 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:29 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:30 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:31 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:32 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:33 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:34 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:35 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:36 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:37 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:38 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:39 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:40 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:41 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:42 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:43 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:44 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:45 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:46 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:47 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:48 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:49 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:50 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:51 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:52 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:53 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:54 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:55 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:56 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:57 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:58 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:59 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:60 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:61 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:62 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:63 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:64 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:65 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:66 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:67 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:68 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:69 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:70 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:71 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:72 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:73 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:74 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:75 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:76 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:77 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:78 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:79 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:80 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:81 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:82 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:83 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:84 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:85 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:86 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:87 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:88 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:89 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:90 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:91 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:92 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:93 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:94 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:95 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:96 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:97 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:98 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:99 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:100 w l \n");
	fprintf(gnuplot, "plot 'DatFile.dat' u 1:101 w l \n");

	pclose(gnuplot);

	system( (char*) "ps2pdf DatFile.ps" );	

	
	delete[] Ham;
	delete[] x_var;
	delete[] eigenvalue;
	delete[] eigenvector;


	return 0;
}