#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>
#include <random>

#include "utils.cpp"

using namespace std;
using namespace arma;

int main(){
    int N = 4;
    
    // mat Lattice = make_sys(N, "random");
    // mat Lattice1 = make_sys(N, "lowest");
    // mat Lattice2 = make_sys(N, "highest");
    // Lattice.print();
    // cout << calc_E(Lattice) << endl;
    // Lattice1.print();
    // cout << calc_E(Lattice1) << endl;
    // Lattice2.print();
    // cout << calc_E(Lattice2) << endl;

    // double e, m, Cv, chi;
    // double T = 1.;

    // mc_cycle(Lattice, T, e, m, Cv, chi);
    // cout << e << endl;
    // mat Lattice3 = make_sys(N, "high");

    // double e;
    // double m;
    // double Cv;
    // double chi;

    // mc_run(N, 600, 100., e, m, Cv, chi, "random");
    // cout << e << " " << m << " " << Cv << " " << chi << endl;

    mat Data = mc_run_culm(N, 10, 1., "lowest");
    Data.print();

    // mat Prob = mc_e_prob(Lattice, 100., 150);
    // Prob.print();

    return 0;
}