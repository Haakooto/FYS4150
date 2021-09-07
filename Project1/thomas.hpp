#ifndef THOMAS_H
#define THOMAS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

void init_arrays(int, double, double *, double *, double *, double *, double *, double *);
void Fwd_Bkwd_sub(int, double *, double *, double *, double *);
void Thomas(int);

#endif