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


void write_cube_to_file(arma::cube C, arma::vec t, string fname, int frame_rate=1){
	/*
	Write position and velocity history of all particles in trap to file.

	Arguments:
		C: arma::cube
			object with positions and velocities
		t: arma::vec
			object with times
		fname: string
			filename. Must include path and extension
		frame_rate: int
			Write only every frame_rate to file. Defualts to write all
	Returns:
		None
	Outputs:
		file with filename
	*/
	ofstream out;
	out.open(fname);
	out << "time particle x y z vx vy vz\n";
	out << fixed << setprecision(8);
	for (int i=0; i < C.n_cols; i++){  // loop over particles
		for (int j=0; j < C.n_slices; j+=frame_rate){  // loop over time
			out << t(j) << " " << i + 1;
			for (int k=0; k < C.n_rows; k++){ // loop over coordinate
				out << " " << C.slice(j).row(k)(i);
			}
			out << endl;
		}
	}
	out.close();
}

void write_errors_to_file(arma::vec err, arma::vec t, string filename){
	/*
	Write relative error in every time-step to file

	Arguments:
		err: arma::vec
			Vector with all relative errors.
			First value is max absolute error. error is 0 here anyways
		t: armma::vec
			times
		filename: string
			filename. Must include path and extension
	Returns:
		None
	Outputs:
		file with filename
	*/
    ofstream out;
	cout << filename << endl;
    out.open(filename);
    out << "t err" << endl;
    out << setprecision(8);
    for (int i=0; i < err.n_rows; i++){
        out << t(i) << " " << err(i) << endl;
    }
out.close();
}


void single_particle_endurace(){  // Ex9p1
	double T = 100;
	double h = 0.005;

	Particle p = Particle(arma::vec({0, 0, 10}), arma::vec({0, 0, 0}), m, q);
	PenningTrap Trap = PenningTrap(b, v, d, false);
	Trap.insert_particles(p);

	Trap.simulate(T, h);

	write_cube_to_file(Trap.get_history(), Trap.get_time(), "outputs/oneP_endurance.txt");
}

void single_particle_errors(string method="RK4"){  // Ex9p5/6
	double x0 = 10;
    double z0 = 10;
    double y_v0 = 10;
    double T_tot = 5;
	int N = 6;

	Particle p = Particle(arma::vec({x0, 0, z0}), arma::vec({0, y_v0, 0}), m, q);
	PenningTrap Trap = PenningTrap(b, v, d, false);
	Trap.insert_particles(p);

	for (int i = 0; i < N; i++){
		arma::vec errs = Trap.analytic_analysis(T_tot, pow(10, -i), method);
		cout << errs(0) << endl;
		write_errors_to_file(errs, Trap.get_time(), "outputs/rel_errors_" + method + "_neglog10dt_" + to_string(i) + ".txt");
	}
}

void two_particle(){  // Ex9p2/3/4
	double T = 100;
	double h = 0.01;

	PenningTrap Trap = PenningTrap(b, v, d, true);
	Particle p1 = Particle(arma::vec({20, 0, 10}), arma::vec(3, arma::fill::randn), m, q);
	Particle p2 = Particle(arma::vec({-20, 0, -10}), arma::vec(3, arma::fill::randn), m, q);

	Trap.insert_particles(p1);
	Trap.insert_particles(p2);
	Trap.simulate(T, h);
	write_cube_to_file(Trap.get_history(), Trap.get_time(), "outputs/twoP_ppi.txt");

	Trap.ppi = false;
	Trap.simulate(T, h);
	write_cube_to_file(Trap.get_history(), Trap.get_time(), "outputs/twoP_no_ppi.txt");
}

void wallpaper(){  // For fun. Run make_cool_wallpaper in plot.py after. Look at the resulting figure from above
	double T = 10000;
	double h = 0.01;

	PenningTrap Trap = PenningTrap(b, v, d, true);
	Particle p1 = Particle(arma::vec({20, 0, 20}), arma::vec(3, arma::fill::randn), m, q);
	Particle p2 = Particle(arma::vec({-20, 0, -20}), arma::vec(3, arma::fill::randn), m, q);

	Trap.insert_particles(p1);
	Trap.insert_particles(p2);
	Trap.simulate(T, h);
	write_cube_to_file(Trap.get_history(), Trap.get_time(), "outputs/twoP_ppi_tall.txt");

	Trap.ppi = false;
	Trap.simulate(T, h);
	write_cube_to_file(Trap.get_history(), Trap.get_time(), "outputs/twoP_no_ppi_tall.txt");
}

double V(double t, double V0, double f, double w){
	// time-dependent potential, eq 21 in project description
	return V0 * (1 + f * cos(w * t));
}

void broad_freq_search(){  // Ex10p1
	double T = 500;
	double timestep = 0.01;

	int N = 100;  // number of particles
	double sd = 0.05;  // factor difference in d
	double sv = 4000;  // factor difference in v0

	vector<double> amps = {0.1, 0.4, 0.7};
	double w_min = 0.2;
	double w_max = 2.51;
	double w_step = 0.02;

	// Calculate some parameters
	double w0 = q * b / m;
	double wz_sq = 2 * q * (v / sv) / m * pow(d * sd, -2);
	double w_p = (w0 + pow(pow(w0, 2) - 2 * wz_sq, 0.5)) / 2;
	double w_m = (w0 - pow(pow(w0, 2) - 2 * wz_sq, 0.5)) / 2;

	ofstream out;
	out.open("outputs/broad_freq_search_new.txt");
	out << " wz_sq w_min w_plus\n";
	out << wz_sq << " " << w_m << " " << w_p << endl;
	out << "ampl wV fracRem\n";
	out << fixed << setprecision(6);

	// Make a PenningTrap with time-dep Efield. Set to dummy func, returning t
	PenningTrap TimeTrap = PenningTrap(b, [](double t){return t;}, d * sd, false);
	TimeTrap.insert_particles(N, m, q); // insert N particles
	for (double f: amps){
		cout << "f = " << f << endl;
		for (double wV = w_min; wV <= w_max; wV += w_step){
			// Set actual Efield func
			TimeTrap.set_tEfield([&](double t){return V(t, v / sv, f, wV);});
			// simulate function resets the particles, so do not have to reinitialize trap, can just restart simulation with new Efield func
			TimeTrap.simulate(T, timestep);

			double fraq = (double)(N - TimeTrap.escaped()) / N;

			// write to file
			out << f << " " << wV << " " << fraq << endl;
			cout << " wV = " << wV << " ratio remaining: " << fraq << endl;
		}
	}
	out.close();
}


void narrow_freq_search(){  // Ex10p2
    double T = 500;
	double timestep = 0.01;

	int N = 100;  // number of particles
	double sd = 0.05;  // factor difference in d
	double sv = 4000;  // factor difference in v0

	vector<double> amps = {0.1};
	double w_min = 0.35;
	double w_max = 0.55;
	double w_step = 0.005;
    double w0 = q * b / m;

	vector<string> ppis = {"no", "with"};

	// Make a PenningTrap with time-dep Efield. Set to dummy func, returning t
	PenningTrap TimeTrap = PenningTrap(b, [](double t){return t;}, d * sd, false);
	TimeTrap.insert_particles(N, m, q); // insert N particles

	ofstream out;
	for (string m: ppis){
		out.open("outputs/narrow_freq_search_" + m + "_ppi.txt");
		out << "ampl wV fracRem\n";
		out << fixed << setprecision(6);

		for (double f: amps){
			cout << "f = " << f << endl;
			for (double wV = w_min; wV <= w_max; wV += w_step){
				// Set actual Efield func
				TimeTrap.set_tEfield([&](double t){return V(t, v / sv, f, wV);});

				// simulate function resets the particles, so do not have to reinitialize trap, can just restart simulation with new Efield func
				TimeTrap.simulate(T, timestep);

				// Calculate some parameters
				double fraq = (double)(N - TimeTrap.escaped()) / N;

				// write to file
				out << f << " " << wV << " " << fraq << endl;
				cout << m << " interactions: " <<  "wV = " << wV << " ratio remaining: " << fraq << endl;
			}
		}
		out.close();
		TimeTrap.ppi = true;
	}
}


void ex10_particle_track(){
    double T = 500;
    double timestep = 0.005;

    double sd = 0.05;  // factor difference in d
    double sv = 4000;  // factor difference in v0

    double f = 0.4;
    double wV = 0.5;

    PenningTrap TimeTrap = PenningTrap(b, [&](double t){return V(t, v / sv, f, wV);}, d * sd, false);
    Particle p = Particle(arma::vec({50, 50, 50}), arma::vec({0, 50, 30}), m, q);
	TimeTrap.insert_particles(p);

    TimeTrap.simulate(T, timestep);

	write_cube_to_file(TimeTrap.get_history(), TimeTrap.get_time(), "outputs/ex10_TimeTrap_particle_track_f0.4_w0.49.txt");

    PenningTrap RegularTrap = PenningTrap(b, v/sv, d * sd, false);
    RegularTrap.insert_particles(p);

    RegularTrap.simulate(T, timestep);

    write_cube_to_file(RegularTrap.get_history(), RegularTrap.get_time(), "outputs/ex10_RegularTrap_particle_track_f0.4_w0.49.txt");
}


void ex10_particles_escape(){
    double T = 500;
    double timestep = 0.005;
    int N = 100;

    double sd = 0.05;  // factor difference in d
    double sv = 4000;  // factor difference in v0

    double f = 0.7;
    double wV = 2.38;

    PenningTrap TimeTrap = PenningTrap(b, [&](double t){return V(t, v / sv, f, wV);}, d * sd, false);
	TimeTrap.insert_particles(N, m, q);

    TimeTrap.simulate(T, timestep);

    double fraq = (double)(N - TimeTrap.escaped()) / N;

    cout << "Particles remaining: " << fraq*100 << endl;


    }




void run_all_experiments(){
	single_particle_endurace();  // Ex9p1
	single_particle_errors();  // Ex9p5/6
	single_particle_errors("Euler");  // Ex9p5/6 Euler
	two_particle();  // Ex9p2/3/4
	wallpaper();
	broad_freq_search(); // Ex10p1
    narrow_freq_search();  // Ex10p2
	ex10_particle_track();
	ex10_particles_escape();
}


int main() {
	run_all_experiments();
	return 0;
}
