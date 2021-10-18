#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <complex>
#include <armadillo>
#include <time.h>

using namespace std;
#include "PenningTrap.hpp"

// Define physical constants
extern const double ke = 1.38935333 * pow(10, 5);  // u (mu m)^3 / (mu s)^2 / e^2

void part1(double amplitude){
    double d = 0.05 * pow(10, 4);
    double V0 = 0.0025 * pow(10, 8);
    double B = 96.5;
    int q = 1; // e
    double m = 40.078; // u Ca+
    int N = 100;
    double f = amplitude;
    double T_tot = 500;
    double timestep = 0.05;
    double w_step = 0.2;
    int steps = int(2.5/w_step);


    arma::vec fraction(steps, arma::fill::zeros);
    arma::vec w_v_vec(steps, arma::fill::zeros);

    arma::vec w_plus(steps, arma::fill::zeros);
    arma::vec w_min(steps, arma::fill::zeros);
    arma::vec w_z_sq(steps, arma::fill::zeros);

    double w_V = w_step;

    for (int i=0; i < steps; i++){
        double w_0 = q * B / m;
        //cout << w_V << endl;
        w_z_sq(i) = 2*q * V0 * (1 + f * cos(w_V * T_tot))/(m * pow(d, 2));
        w_plus(i) = (w_0 + pow(pow(w_0, 2) - 2 * w_z_sq(i), 0.5))/2;
        w_min(i) = (w_0 - pow(pow(w_0, 2) - 2 * w_z_sq(i), 0.5))/2;

        w_v_vec(i) = w_V;

        PenningTrap Trap = PenningTrap(B, V0, d, false, f, w_V);
        Trap.insert_particles(N, m, q);
        Trap.simulate(T_tot, timestep, "RK4");

        int escaped = Trap.escaped();

        fraction(i) = (double)(N - escaped)/N;

        w_V += w_step;
    }


    ofstream out;
    string filename = "outputs/time_dep_V_amp_" + to_string((int)(10*f)) + ".txt";
    out.open(filename);
    out << fixed << setprecision(6);
    out << "w_V frac w_z_sq w_min w_plus";
    out << endl;

    for (int i=0; i < steps; i++){
        out << w_v_vec(i) << " " << fraction(i) << " " << w_z_sq(i) << " " << w_min(i) << " " << w_plus(i);
        out << endl;
    }
    out.close();
    cout << "Done writing" << endl;
}



void run_all_experiment(){
    double amplitude = 0.1;
    part1(amplitude);
}




int main(){

    double amplitude = 0.7;
    part1(amplitude);

    return 0;
}
