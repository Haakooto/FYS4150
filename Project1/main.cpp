#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double f(double x){return 1 - (1 - exp(-10)) * x - exp(-10 * x);}

int main(){
    double x0 = 0; 
    double x1 = 1;
    int N = 1001;
    double h = (x1 - x0) / N;
    double x[N];
    double u[N];

    for (int i = 0; i < N; i++){
        x[i] = i * h;
        u[i] = f(x[i]);
    }

    ofstream out;
    out.open("data.txt");
    out << "x u\n";
    out << fixed;
    for (int i = 0; i < N; i++){
        out << setprecision(8) << x[i] << " " << u[i] << "\n";
    }
    out.close();
    return 0;
}
