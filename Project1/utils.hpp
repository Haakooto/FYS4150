#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x){return 100 * exp(-10 * x);}
double analytic_sol(double x){return 1 - (1 - exp(-10)) * x - exp(-10 * x);}

void write_to_file(int n, double *x, double *u, double *g, double *aerr, double *rerr){
    ofstream out;
    out.open("datas/data_" + to_string(n) + ".txt");
    out << "x u v abs_err rel_err\n";
    out << fixed;
    for (int i = 0; i < pow(10, n) + 1; i++){
        out << setprecision(10) << x[i];
        out << " " << u[i];
        out << " " << g[i];
        out << " " << aerr[i];
        out << " " << rerr[i];
        out << endl;
    }
    out.close();
}

double *rel_err(int N, double*u, double *v){
    double *error = new double[N];
    for (int i = 0; i < N; i++){
        error[i] = log10((u[i] - v[i]) / u[i]);
    }
    return error;
}

double *abs_err(int N, double*u, double *v){
    double *error = new double[N];
    for (int i = 0; i < N; i++){
        error[i] = log10(abs(u[i] - v[i]));
    }
    return error;
}

#endif