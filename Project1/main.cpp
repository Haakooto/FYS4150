#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x){return 100 * exp(-10 * x);}
double analytic_sol(double x){return 1 - (1 - exp(-10)) * x - exp(-10 * x);}

void write_to_file(int n, double *x, double *u, double *g){
    ofstream out;
    out.open("data_" + to_string(n) + ".txt");
    out << "x u v\n";
    out << fixed;
    for (int i = 0; i < pow(10, n) + 1; i++){
        out << setprecision(10) << x[i] << " " << u[i] << " " << g[i]  << "\n";
    }
    out.close();
}

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

double *RelErr(int N, double*u, double *g){
    double *error = new double
    [N];
    for (int i = 0; i < N; i++){
        error[i] = log10((u[i] - g[i]) / u[i]);
    }
    return error;
}

double *AbsErr(int N, double*u, double *g){
    double *error = new double[N];
    for (int i = 0; i < N; i++){
        error[i] = log10(abs(u[i] - g[i]));
    }
    return error;
}

void Thomas(int n){
    int N = pow(10, n) + 1;
    double x[N];
    double u[N];
    double a[N];
    double b[N];
    double c[N];
    double g[N];
    
    double x0 = 0;
    double x1 = 1;
    double h = (x1 - x0) / (N - 1);

    init_arrays(N, h, x, u, a, b, c, g);
    Fwd_Bkwd_sub(N, a, b, c, g);

    double *rel_err = RelErr(N, u, g);
    double *abs_err = AbsErr(N, u, g);
    
    write_to_file(n, x, u, g);

    delete[] x;
    delete[] u;
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] g;

}


int main(int argc, char *argv[]){
    double x0 = 0;
    double x1 = 1;
    int n;

    n = atoi(argv[1]);
    Thomas(n);

    return 0;
}
