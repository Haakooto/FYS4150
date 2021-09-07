#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

long double f(long double x){return 100 * exp(-10 * x);}
long double analytic_sol(long double x){return 1 - (1 - exp(-10)) * x - exp(-10 * x);}

void write_to_file(int n, long double *x, long double *u, long double *g, long double *aerr, long double *rerr){
    ofstream out;
    out.open("datas/data_" + to_string(n) + ".txt");
    out << "x u v abs_err rel_err\n";
    out << fixed;
    for (int i = 0; i < pow(10, n) + 1; i++){
        out << setprecision(14) << x[i];
        out << " " << u[i];
        out << " " << g[i];
        out << " " << aerr[i];
        out << " " << rerr[i];
        out << endl;
    }
    out.close();
}

long double *rel_err(int N, long double*u, long double *v){
    long double *error = new long double[N];
    for (int i = 0; i < N; i++){
        error[i] = log10((u[i] - v[i]) / u[i]);
    }
    return error;
}

long double *abs_err(int N, long double*u, long double *v){
    long double *error = new long double[N];
    for (int i = 0; i < N; i++){
        error[i] = log10(abs(u[i] - v[i]));
    }
    return error;
}

#endif