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
// Define Ca+ ion properties
const int q = 1; // e
const double m = 40.078; // u Ca+
// Define standard trap properties
const double b = 96.5;  // u / (mu s) / e
const double v = 9.65 * pow(10, 8);  // u (mu m)^2 / (mu s)^2 / e
const double d = pow(10, 4); // mu m


class TimePotential{
  public:
	double V0, f, wV;
	TimePotential(double strength, double ampl, double freq){
		V0 = strength;
		f = ampl;
		wV = freq;
	}
	double call(double t){
		return V0 * (1 + f * cos(wV * t));
	}
};


void write_cube_to_file(arma::cube C, arma::vec t, string fname, int frame_rate=1){
	ofstream out;
	out.open(fname);
	out << "time particle x y z\n";
	out << fixed << setprecision(8);
	for (int i=0; i < C.n_slices; i++){  // loop over particles
		for (int j=0; j < C.n_rows; j+=frame_rate){  // loop over time
			out << t(j) << " " << i + 1;
			for (int k=0; k < C.n_cols; k++){ // loop over coordinate
				out << " " << C.slice(i).row(j)(k);
			}
			out << endl;
		}
	}
	out.close();
}


void write_analytic_solution_to_file(arma::mat R, arma::vec t, string filename){
    ofstream out;
    out.open(filename);
    out << "t x y z" << endl;
    out << fixed << setprecision(8);
    for (int i=0; i < R.n_rows; i++){
        out << t(i);
        out << " " << R(i, 0) << " " << R(i, 1) << " " << R(i, 2);
        out << endl;
    }
out.close();
}

double f(double t){
	return sin(t);
}


void z_movement(){
	double T = 100;
	double h = 0.005;

	Particle p = Particle(arma::vec({0, 0, 10}), arma::vec({0,0,0}), m, q);
	PenningTrap Trap = PenningTrap(b, v, d, false);
	Trap.insert_particles(p);

	Trap.simulate(T, h);

	write_cube_to_file(Trap.get_history(), Trap.get_time(), "outputs/z_movement.txt");
}

void ppi_comparison(){
	double T = 10;
	double h = 0.005;

	Particle p1 = Particle(arma::vec(3, arma::fill::randu), arma::vec(3, arma::fill::randn), m, q);
	Particle p2 = Particle(arma::vec(3, arma::fill::randu), arma::vec(3, arma::fill::randn), m, q);
	
	PenningTrap Trap = PenningTrap(b, v, d, true);
	PenningTrap Trap_no_ppi = PenningTrap(b, v, d, false);

	Trap.insert_particles(p1);
	Trap.insert_particles(p2);
	Trap.simulate(T, h);

	Trap_no_ppi.insert_particles(p1);
	Trap_no_ppi.insert_particles(p2);
	Trap_no_ppi.simulate(T, h);

	write_cube_to_file(Trap.get_history(), Trap.get_time(), "outputs/ppi_comparison.txt");
	write_cube_to_file(Trap_no_ppi.get_history(), Trap_no_ppi.get_time(), "outputs/ppi_comparison_no_ppi.txt");
}



void single_particle(){
	double x0 = 10;
    double z0 = 10;
    double y_v0 = 10;
    double T_tot = 1;
    double timestep = 0.00005;

	Particle p = Particle(arma::vec({x0,0,z0}), arma::vec({0,y_v0,0}), m, q);
	PenningTrap Trap = PenningTrap(b, v, d, false);
	Trap.insert_particles(p);

	Trap.simulate(T_tot, timestep);
	Trap.analytic(T_tot, timestep, x0, z0, y_v0);


}


void experiments(){
	z_movement();  // first point in P9
	ppi_comparison();  // second point in P9
}



int main() {
	// define particle and trap properties.


	// TimePotential TP = TimePotential(v, 1, 1); // not sure what f and wV should be

    // PenningTrap TimeTrap = PenningTrap(b, (*f), d, true);
    // PenningTrap TimeTrap = PenningTrap(b, (*TP.call), d, true);

	//Particle p2 = Particle(arma::vec({0,0,-1}), arma::vec({0,0,0}), 1, 1);
	// Particle p3 = Particle(arma::vec({0,0,-200}), arma::vec({0,0,0}), 1, 1);
	// Particle p4 = Particle(arma::vec({0,0,200}), arma::vec({0,0,0}), 1, 1);

	//Trap.insert_particles(100, m, q);
	Trap.insert_particles(p1);
	// Trap.insert_particles(p2);
	// Trap.insert_particles(p3);
	// Trap.insert_particles(p4);

	//clock_t t1 = clock();
	Trap.simulate(T_tot, timestep, "Euler");
	Trap.analytic(T_tot, timestep, x0, z0, y_v0);
    //clock_t t2 = clock();
    //double time = ((double)(t2 - t1) / CLOCKS_PER_SEC);
	//cout << time << endl;
	write_cube_to_file(Trap.get_history(), Trap.get_time(), "test.txt");
	write_analytic_solution_to_file(Trap.get_asol(), Trap.get_time(), "analytic_solution.txt");
	return 0;
}
