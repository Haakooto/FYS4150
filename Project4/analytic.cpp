#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>
#include <random>

#include "ising_model.cpp"

using namespace std;
using namespace arma;


int main(int argc, char* argv[])
/*
Calculates the analytical solutions for L=2, and the average of R numerical solutions over M Monte Carlo cycles.
The results are printed in the terminal, either in text form or as a Latex table.

Arguments:
    T: double
        Temperature
    M: int
        No. of Monte Carlo cycles
    R: int
        No. of runs to take the average over
*/
{
    int M, R;
    double T;
    bool ugly_out = false;
    if (argc != 4 && argc != 5){
        cout << "Bad usage! This program takes three parameters: ";
        cout << "temperature, number of monte carlo cycles, and rounds to average over \n";
        return 1;
    } else {
        T = atof(argv[1]);
        M = atoi(argv[2]);
        R = atoi(argv[3]);
    }
    if (argc == 5){
        ugly_out = true;
    }


    int L = 2;
    int N = L*L;
    string method = "random";
    int burnin = 0;

    // vector with calculated quantities. See docstring of multi_mc for content
    arma::vec data(8, arma::fill::zeros);

    multi_mc(L, M, R, T, data, method, burnin, true);

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

    if (ugly_out) {  // If called from python, just print outputs
        arma::rowvec out(12);
        out(0) = e_a; out(1) = m_a; out(2) = Cv_a; out(3) = chi_a;
        out.cols(4, 11) = data.t();
        out.raw_print(std::cout);

    } else {  // if not called from python, make nice prints
        cout << "" << endl;
        cout << "Analytic solutions" << endl;
        cout << "Average energy per spin: " << e_a  << endl;
        cout << "Average magnetisation per spin: " << m_a << endl;
        cout << "Specific heat capacity: " << Cv_a << endl;
        cout << "Susceptibility: " << chi_a << endl;
        cout << "" << endl;

        cout << "" << endl;
        cout << "Numerical solutions" << endl;
        cout << "Average energy per spin: " << data(0) << " ± " << data(4) << endl;
        cout << "Average magnetisation per spin: " << data(1) << " ± " << data(5) << endl;
        cout << "Specific heat capacity: " << data(2) << " ± " << data(6) << endl;
        cout << "Susceptibility: " << data(3) << " ± " << data(7) << endl;
    }

    return 0;
}
