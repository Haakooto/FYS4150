#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <armadillo>
#include <time.h>

using namespace std;

const double ke = 1.38935333 * pow(10, 5);

class Particle{
  public:
	double m;
	int q;
	arma::vec r, v;

	Particle(arma::vec position, arma::vec velocity, double mass, int charge){
		r = position;
		v = velocity;
		m = mass;
		q = charge;
	}

	void print(){
		cout << "Particle at ( ";
		for (int i=0; i<3; i++){cout << r(i) << " ";}
		cout << ") with velocity ( ";
		for (int i=0; i<3; i++){cout << v(i) << " ";}
		cout << ") with mass " << m << " and charge " << q << endl;
	}

	// void advance(arma::vec dr, arma::vec dv){
	// 	r += dr;
	// 	v += dv;
	// }
};

class PenningTrap{
  public:
	int N = 0;
	double B0, V0, d;
	double (*tV0)(double);
	bool ppi, time_dep_V;
	vector<Particle> particles;
	arma::cube r;
	arma::mat v;

	PenningTrap(double Bfield, double Efield, double length, bool particle_particle=false){
		B0 = Bfield;
		V0 = Efield;
		d = length;
		time_dep_V = false;
		ppi = particle_particle;
	}

	PenningTrap(double Bfield, double(*tEfield)(double), double length, bool particle_particle=true){
		B0 = Bfield;
		tV0 = tEfield;
		d = length;
		time_dep_V = true;
		ppi = particle_particle;
	}

	void insert_particles(vector<Particle> P){
		for (Particle p: P){
			particles.push_back(p);
			N++;
		}
	}

	void insert_particles(Particle p){
		particles.push_back(p);
		N++;
	}

	void insert_particles(int n, double q, double m){ // inserting n identical particles uniformly in sphere
		// https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
		for (int i=0; i<n; i++){
			arma::vec r(3, arma::fill::randn);
			arma::vec u(1, arma::fill::randu);
			r *= d * pow(u(0), 3) / arma::norm(r);
			particles.push_back(Particle(r, arma::zeros(3), m, q));
			N++;
		}
	}

	double get_Efield_at_time(double t){
		if (time_dep_V) {return tV0(t);}
		else {return V0;}
	}

	arma::mat sum_particles_forces(arma::mat ri){
		/*
		Calculates all particle-particle forces
		ri has shape 3, N
		*/
		arma::mat Eforce(N, 3, arma::fill::zeros);
		if (!ppi) {return Eforce.t();}   //ignore particle-particle interactions

		else {
            for (int i=0; i < N; i++){
				for (int j=0; j < i; j++){
                    arma::vec diff_ri = ri.col(i)-ri.col(j);
                    double norm = arma::norm(diff_ri);
                    arma::vec F = ke * Q(i) * Q(j) * diff_ri * pow(norm, -3);
					Eforce.row(i) += F.t() / M(i);
					Eforce.row(j) -= F.t() / M(j);
				}
			}
		} return Eforce.t();  // 3 x N
	}

	double dt;
	int nT;
	arma::vec t;  // time-vector
    arma::rowvec Q; //list of all charges	
    arma::rowvec M; //list of all masses

	void simulate(double T, double timestep){
		dt = timestep;
		nT = (int)(T / dt) + 1;
		t = arma::vec(nT, arma::fill::zeros);
		r = arma::cube(nT, 3, N);  // nT timesteps x 3dim x N particles
		v = arma::mat(3, N);  // 3dim x N particles
		
		Q = arma::rowvec(N);
		M = arma::rowvec(N);

		// initiate simulation
		for (int i=0; i < N; i++){
			Particle p = particles[i];
			r.slice(i).row(0) = p.r.t();  // particle i at time 0
			v.col(i) = p.v;  
            Q(i) = p.q;
            M(i) = p.m;
		}
		// start simulation
		for (int i=0; i < nT - 1; i++){
			RK4(i, t(i));
			t(i + 1) = (i + 1) * dt;
		}
	}

	void RK4(int i, double t){
		double h = dt / 2;
		arma::cube u, k1, k2, k3, k4;
		u = arma::cube(3, N, 2);  // 3dim x N particles x 2 phases (position and velocity)

		// wacky bugfix.
		arma::mat r_(r.row(i));
		if (N == 1){ u.slice(0) = r_.t();}  // position at time i
		else { u.slice(0) = r_;}  // position at time i		
		u.slice(1) = v;  // velocity

		k1 = advance(t, u);
		k2 = advance(t + h, u + h * k1);
		k3 = advance(t + h, u + h * k2);
		k4 = advance(t + 2 * h, u + h * k3);
		u += h * (k1 + 2 * k2 + 2 * k3 + k4) / 3;

		r.row(i + 1) = u.slice(0);
		v = u.slice(1);

		// check for escapees
		for (int p=0; p < N; p++){
			if (arma::norm(r.slice(p).row(i + 1)) > d){
				Q(p) = 0;
			}
		}

	}

	arma::cube advance(double time, arma::cube u){
		arma::cube du(size(u));  // change in u
		arma::mat F = sum_particles_forces(u.slice(0));  // calculate ppi from poistions
		du.slice(0) = u.slice(1);  // set change in pos to vel

		Efield(u.slice(0), time);
		Bfield(u.slice(1));

		//update velocities for x, y, z for all particles
		du.slice(1).row(0) = u.slice(0).row(0) + u.slice(1).row(1) + F.row(0);  
		du.slice(1).row(1) = u.slice(0).row(1) + u.slice(1).row(0) + F.row(1);
		du.slice(1).row(2) = u.slice(0).row(2) + F.row(2);

		return du;  // 3dim x N particles x 2 phases (position and velocity)
	}

	void Efield(arma::mat &pos, double t){
		double V = get_Efield_at_time(t);
		pos.row(0) %= Q / M * pow(d, -2) * V;
		pos.row(1) %= Q / M * pow(d, -2) * V;
		pos.row(2) %= -2 * Q / M * pow(d, -2) * V;
	}

	void Bfield(arma::mat &vel){
		vel.row(0) %= Q * B0 / M;
		vel.row(1) %= -Q * B0 / M;
	}
};


double f(double t){
	return sin(t);
}

void write_cube_to_file(arma::cube C, arma::vec t, string fname){
	ofstream out;
	out.open(fname);
	out << "t";
	for (int i=0; i < C.n_slices; i++){
		out << " x" + to_string(i) + " y" + to_string(i) + " z" + to_string(i);
	}
	out << endl;
	out << fixed << setprecision(8);
	for (int i=0; i < C.n_rows; i++){
		out << t(i);
		for (int j=0; j < C.n_slices; j++){
			for (int k=0; k < C.n_cols; k++){
				out << " " << C.slice(j).col(k)(i);  
			}
		}
		out << endl;
	}
	out.close();
}


int main() {
	Particle p1 = Particle(arma::vec({0,0,1}), arma::vec({0,0,0}), 1, 1);
	Particle p2 = Particle(arma::vec({0,0,-1}), arma::vec({0,0,0}), 1, 1);
	// Particle p3 = Particle(arma::vec({0,0,-200}), arma::vec({0,0,0}), 1, 1);
	// Particle p4 = Particle(arma::vec({0,0,200}), arma::vec({0,0,0}), 1, 1);
	double b = 9.65;
	double v = 9.65 * pow(10, 8);
	double d = pow(10, 4);

    PenningTrap P = PenningTrap(b, v, d, true);
	P.insert_particles(100, 1, 1);
	// P.insert_particles(p1);
	// P.insert_particles(p2);
	// P.insert_particles(p3);
	// P.insert_particles(p4);
	// P.r.print();
	clock_t t1 = clock();
	P.simulate(1, 0.0005);
    clock_t t2 = clock();
    double time = ((double)(t2 - t1) / CLOCKS_PER_SEC);
	cout << time << endl;
	// P.r.print();
	// P.r.reshape(0,2,1).save("r.csv", arma::file_type::arma_ascii);
	write_cube_to_file(P.r, P.t, "test.txt");
	// arma::cube pos(P.r);
	cout << P.N - arma::sum(P.Q) << endl;
	// pos.reshape(P.r.n_rows,P.r.n_slices,P.r.n_cols);
	// pos.save("r.csv", arma::file_type::arma_ascii);
	// arma::mat a(4, 4, arma::fill::randn);
	// // arma::vec b(4, arma::fill::randn);
	// a.print();
	// cout << arma::norm(a) << endl;
	// b.print();
	// (a * b).print();

	// // cout << P.get_Efield_at_time(0.2) << endl;
	// cout << P.sum_particles_forces() << endl;
	// arma::mat U0(3, 2, arma::fill::randn);
	// U0.print();
	// arma::mat U(size(U0));
	// U.print();
	// // arma::cube R(2, 3, 4);
	// arma::cube U(4, U0.n_rows, U0.n_cols); //, arma::fill:zeros);
	// U.randn();
	// U.print();
	// U.slice(0).row(0).print();
	// bool a = true;
	// if (!a) { cout << "true\n";}
	// else {cout << a << "false\n";}


	return 0;
}
