#ifndef THOMAS_H
#define THOMAS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "utils.hpp"

using namespace std;

void init_arrays(int N, double h, double *x, double *u, double *a, double *b, double *c, double *g){
    for (int i = 0; i < N ; i++){
        x[i] = i * h;
        u[i] = analytic_sol(x[i]);

        a[i] = -1;
        b[i] = 2;
        c[i] = -1;
        g[i] = pow(h, 2) * f(x[i]);
    }
    // Set boundary conditions
    g[0] = 0;
    g[N - 1] = 0;
}

void Fwd_Bkwd_sub(int N, double *a, double *b, double *c, double *g){
    for (int i = 2; i < N - 1; i++){
        b[i] = b[i] - a[i] / b[i - 1] * c[i];
        g[i] = g[i] - a[i] / b[i - 1] * g[i - 1];
    }
    for (int i = N - 2; i > 0; i--){
        g[i] = (g[i] - c[i] * g[i + 1]) / b[i];
    }
}

void Thomas(int n){
    int N = pow(10, n) + 1;

    double *x, *u, *a, *b, *c, *g;
    x = new double[N];
    u = new double[N];
    a = new double[N];
    b = new double[N];
    c = new double[N];
    g = new double[N];

    double x0 = 0;
    double x1 = 1;
    double h = (x1 - x0) / (N - 1);

    cout << "Initializing arrays\n";
    init_arrays(N, h, x, u, a, b, c, g);
    cout << "Performing forward and backward substituion\n";
    Fwd_Bkwd_sub(N, a, b, c, g);

    cout << "Calculating absolute error\n";
    double *aerr = abs_err(N, u, g);
    cout << "Calculating relative error\n";
    double *rerr = rel_err(N, u, g);

    cout << "Writing to file\n";
    write_to_file(n, x, u, g, aerr, rerr);

    delete[] a;
    delete[] x;
    delete[] u;
    delete[] b;
    delete[] c;
    delete[] g;
    delete[] rerr;
    delete[] aerr;
}

#endif