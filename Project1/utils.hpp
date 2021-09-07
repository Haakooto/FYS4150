#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double);
double analytic_sol(double);
void write_to_file(int, double *, double *, double *, double *, double *);
void *abs_err(int, double *, double *);
void *rel_err(int, double *, double *);

#endif