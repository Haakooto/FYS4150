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
    mat Lattice = make_sys(2, "lowest");
    // mat Lattice1 = make_sys(2, "random");
    // mat Lattice2 = make_sys(2, "highest");
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
    // mat Lattice3 = make_sys(2, "high");

    // double e;
    // double m;
    // double Cv;
    // double chi;

    // mc_run_single(2, 600, 100., e, m, Cv, chi, "random");
    // cout << e << " " << m << " " << Cv << " " << chi << endl;

    mat Data = mc_run_culm(2, 10, 100., "random");
    Data.print();

    return 0;
}