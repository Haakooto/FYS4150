#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>
#include <random>

#include "utils.cpp"

using namespace std;
using namespace arma;


int main(int argc, char* argv[])
// Terminal input to run: ./analytic.out  T   M = #no. MC cycles  R = runs to average over
{
    int M, R;
    double T;
    if (argc != 4){
        cout << "Bad usage! This program takes three params: ";
        cout << "temperature, number of monte carlo cycles, and rounds to average over \n";
        return 1;
    } else {
        T = atof(argv[1]);
        M = atoi(argv[2]);
        R = atoi(argv[3]);
    }



    // Test the program against the analytical solution
    int L = 2;
    int N = L*L;
    string method = "random";
    int burnin = 100;

    double E, Mag, CV, CHI;
    double e, m, Cv, chi;


    multi_mc(L, M, R, T, e, m, Cv, chi, method, burnin);



    // Analytic solutions
    double E_a = (-8*sinh(8/T))/(3+cosh(8/T));
    double e_a = E_a/N;
    double E_sq = (64*sinh(8/T))/(3+cosh(8/T));
    double e_sq = E_sq/N;
    double Cv_a = 1/(T*T)*(e_sq - pow(e_a, 2)*N);
    double M_a = (4 + 2*exp(8/T))/(3+cosh(8/T));
    double m_a = M_a/N;
    double M_sq = (8 + 8*exp(8/T))/(3+cosh(8/T));
    double m_sq = M_sq/N;
    double chi_a = (1/T)*(m_sq - pow(m_a, 2)*N);


    cout << "" << endl;
    cout << "Analytic solutions" << endl;
    cout << "Average energy per spin: " << e_a << endl;
    cout << "Average magnetisation per spin: " << m_a << endl;
    cout << "Specific heat capacity: " << Cv_a << endl;
    cout << "Susceptibility: " << chi_a << endl;
    cout << "" << endl;

    cout << "" << endl;
    cout << "Numerical solutions" << endl;
    cout << "Average energy per spin: " << e << endl;
    cout << "Average magnetisation per spin: " << m << endl;
    cout << "Specific heat capacity: " << Cv << endl;
    cout << "Susceptibility: " << chi << endl;



    return 0;
}
