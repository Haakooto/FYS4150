#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <armadillo>
#include <time.h>


#include "Particle.hpp"
#include "PenningTrap.hpp"

using namespace std;

// Define physical constants
extern const double ke = 1.38935333 * pow(10, 5);  // u (mu m)^3 / (mu s)^2 / e^2


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
	for (int i=0; i < C.n_slices; i++){
		for (int j=0; j < C.n_rows; j+=frame_rate){
			out << t(j) << " " << i + 1;
			for (int k=0; k < C.n_cols; k++){
				out << " " << C.slice(i).row(j)(k);
			}
			out << endl;
		}
	}
	out.close();
}

double f(double t){
	return sin(t);
}


int main() {
	// define particle and trap properties.
	int q = 1; // e
	double m = 40.078; // u Ca+
	double b = 96.5;  // u / (mu s) / e
	double v = 9.65 * pow(10, 8);  // u (mu m)^2 / (mu s)^2 / e
	double d = pow(10, 4); // mu m

	// TimePotential TP = TimePotential(v, 1, 1); // not sure what f and wV should be

	// PenningTrap Trap = PenningTrap(b, v, d, true);
    // PenningTrap TimeTrap = PenningTrap(b, (*f), d, true);

	Particle p1 = Particle(arma::vec({0,0,1}), arma::vec({0,0,1}), 1, 1);
	//Particle p2 = Particle(arma::vec({0,0,-1}), arma::vec({0,0,0}), 1, 1);
	// Particle p3 = Particle(arma::vec({0,0,-200}), arma::vec({0,0,0}), 1, 1);
	// Particle p4 = Particle(arma::vec({0,0,200}), arma::vec({0,0,0}), 1, 1);

	//Trap.insert_particles(100, 1, 1);
	// Trap.insert_particles(p1);
	// Trap.insert_particles(p2);
	// Trap.insert_particles(p3);
	// Trap.insert_particles(p4);

	//clock_t t1 = clock();
	// Trap.simulate(1, 0.0005);
    //clock_t t2 = clock();
    //double time = ((double)(t2 - t1) / CLOCKS_PER_SEC);
	//cout << time << endl;
	// write_cube_to_file(Trap.get_history(), Trap.get_time(), "test.txt");

	return 0;
}
