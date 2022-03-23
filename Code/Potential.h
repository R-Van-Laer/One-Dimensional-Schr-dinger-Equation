#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;


double QHO(double x, double *par);

double Tunn(double x, double *par);

double Well(double x, double *par);

double Free(double x, double *par);

double Const(double x, double *par);

double Step(double x, double *par);


#endif
