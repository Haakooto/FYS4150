#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

long double f(long double x){return 100 * exp(-10 * x);}
long double analytic_sol(long double x){return 1 - (1 - exp(-10)) * x - exp(-10 * x);}

void write_full(int n, long double *x, long double *u, long double *g, long double *aerr, long double *rerr, string s = "_"){
    ofstream out;
    out.open("datas/full" + s + to_string(n) + ".txt");
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

void write_limited(int I, int *n, double *time, double *max_rel_err, string s="_", int rep){
    ofstream out;
    out.open("datas/limited" + s + to_string(rep) + ".txt");
    out << "n time max_error\n";
    out << fixed;
    for (int i = 0; i < I; i++){
        out << n[i] << " ";
        out << setprecision(8) << time[i];
        out << " " << max_rel_err[i];
        out << endl;
    }
    out.close();
}

double rel_err(int N, long double *error, long double*u, long double *v){
    double max = 0.;
    for (int i = 0; i < N; i++){
        error[i] = log10(abs((u[i] - v[i]) / u[i]));
        if (error[i] < max) {
            max = error[i];
        }
    }
    return max;
}

void abs_err(int N, long double *error, long double*u, long double *v){
    for (int i = 0; i < N; i++){
        error[i] = log10(abs(u[i] - v[i]));
    }
}

#endif
