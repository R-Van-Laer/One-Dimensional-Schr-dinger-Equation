#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;


double QHO(double x, double *par){

	double k = par[0] * par[1];	
	double y = k * x*x / 2.0;

	return y;
}

double Tunn(double x, double *par){

	double y;

	if( x > par[0] && x < par[1] ){
		y = par[2];
	}
	else{
		y = 0.0;
	}

	return y;
}

double Well(double x, double *par){
	
	double y;

	if( x > par[0] && x < par[1] ){
		y = 0.0;
	}
	else{
		y = par[2];
	}

	return y;
}

double Free(double x, double *par){

	double y = 0.0;

	return y;
}

double Const(double x, double *par){

	double y = par[0];

	return y;
}

double Step(double x, double *par){

	double y;

	if(x < par[0]){
		y = 0.0;
	}
	else{
		y = par[1];
	}

	return y;
}