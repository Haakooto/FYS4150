#ifndef THOMAS_H
#define THOMAS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <time.h>
#include "utils.hpp"

using namespace std;

void init_arrays(int N, long double h, long double *x, long double *u, long double *a, long double *b, long double *c, long double *g){
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

void Fwd_Bkwd_sub(int N, long double *a, long double *b, long double *c, long double *g){
    for (int i = 2; i < N - 1; i++){
        b[i] = b[i] - a[i] / b[i - 1] * c[i - 1];
        g[i] = g[i] - a[i] / b[i - 1] * g[i - 1];
    }
    for (int i = N - 2; i > 0; i--){
        g[i] = (g[i] - c[i] * g[i + 1]) / b[i];
    }
}

void Thomas(int n, int cap, double *duration, double *max_rel_err){
    int N = pow(10, n) + 1;

    long double *x, *u, *a, *b, *c, *g, *aerr, *rerr;
    x = new long double[N];
    u = new long double[N];
    a = new long double[N];
    b = new long double[N];
    c = new long double[N];
    g = new long double[N];
    aerr = new long double[N];
    rerr = new long double[N];

    long double x0 = 0;
    long double x1 = 1;
    long double h = (x1 - x0) / (N - 1);

    init_arrays(N, h, x, u, a, b, c, g);

    clock_t t1 = clock();
    Fwd_Bkwd_sub(N, a, b, c, g);
    clock_t t2 = clock();
    duration[n - 1] = ((double)(t2 - t1) / CLOCKS_PER_SEC);

    abs_err(N, aerr, u, g);
    max_rel_err[n - 1] = rel_err(N, rerr, u, g);

    if (n < cap){
        cout << "Writing to file\n";
        write_full(n, x, u, g, aerr, rerr);
    }

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