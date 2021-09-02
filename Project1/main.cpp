#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x){return 100 * exp(-10 * x);}
double analytic_sol(double x){return 1 - (1 - exp(-10)) * x - exp(-10 * x);}

// void make_analytical_data(int N; double h; **double x; **double analytical_u) {
//     for (int i = 0; i < N; i++){
//         x[i] = i * h;
//         analytical_u[i] = f(x[i]);
//     }
// }

int main(){
    double x0 = 0;
    double x1 = 1;
    int N = 1001;
    double h = (x1 - x0) / (N - 1);
    double x[N];
    double u[N];
    double a[N];
    double b[N];
    double c[N];
    double g[N];

    for (int i = 0; i < N; i++){
        x[i] = i * h;
        u[i] = analytic_sol(x[i]);

        a[i] = -1;
        c[i] = -1;
        b[i] = 2;
        g[i] = pow(h, 2) * f(x[i]);
    }
    // Forward sub
    for (int i = 2; i < N - 1; i++){
        b[i] = b[i] - a[i] / b[i - 1] * c[i];
        g[i] = g[i] - a[i] / b[i - 1] * g[i - 1];
    }
    // for (int i = 0; i < N; i++){
    // cout << g[i] << " ";
    // // cout << c
    // }
    // Backward sub
    g[N-2] = g[N-2] / b[N-2];
    for (int i = N - 3; i > 0; i--){
        g[i] = (g[i] - c[i] * g[i + 1]) / b[i];
    }
    g[N - 1] = 0;
    g[0] = 0;
    ofstream out;
    out.open("data.txt");
    out << "x u v\n";
    out << fixed;
    for (int i = 0; i < N; i++){
        out << setprecision(10) << x[i] << " " << u[i] << " " << g[i]  << "\n";
    }
    out.close();
    return 0;
}
