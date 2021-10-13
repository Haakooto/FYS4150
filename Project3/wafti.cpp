#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>

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

	void advance(arma::vec dr, arma::vec dv){
		r += dr;
		v += dv;
	}

	arma::vec force_from(Particle other){
		arma::vec rirj = r - other.r;
		double norm = arma::norm(rirj);
		return other.q * rirj * pow(norm, -3);
	}
};

class PenningTrap{
  public:
	int N = 0;
	double B0, V0, dsq;
	double (*tV0)(double);
	bool ppi, time_dep_V;
	vector<Particle> particles;
	arma::cube r;
	arma::mat v;

	PenningTrap(double Bfield, double Efield, double length, bool particle_particle=false){
		B0 = Bfield;
		V0 = Efield;
		dsq = pow(length, 2);
		time_dep_V = false;
		ppi = particle_particle;
	}

	PenningTrap(double Bfield, double(*tEfield)(double), double length, bool particle_particle=true){
		B0 = Bfield;
		tV0 = tEfield;
		dsq = pow(length, 2);
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

	void insert_identical_resting_particles(int n, double q, double m){
		// https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
		for (int i=0; i<n; i++){
			arma::vec r(3, arma::fill::randn);
			arma::vec u(1, arma::fill::randu);
			r *= dsq * pow(u(0), 3) / arma::norm(r);  // must use d, not dsq
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
		if (!ppi) {return Eforce;}
		else {
			for (int i=0; i < N; i++){
				Particle alice = particles[i];
				for (int j=0; j < i; j++){
					Particle bob = particles[j];
					arma::vec F = -ke * alice.force_from(bob);
					Eforce.row(i) += F.t() * alice.q / alice.m;
					Eforce.row(j) -= F.t() * alice.q / bob.m;
				}
			}
		} return Eforce;
	}

	double dt;
	int nT;
	arma::vec t;
	arma::vec w0;
	arma::vec wzsq;

	void simulate(double T, double timestep){
		dt = timestep;
		nT = (int)(T / dt);
		t = arma::vec(nT, arma::fill::zeros);
		r = arma::cube(nT, 3, N);  // nT timesteps x 3dim x N particles
		v = arma::mat(3, N);  // 3dim x N particles

		for (int i=0; i < N; i++){
			Particle p = particles[i];
			r.slice(i).row(0) = p.r.t();
			v.col(i) = p.v;
			w0(i) = p.q * B0 / p.m;
			wzsq(i) = 2 * p.q / p.m / dsq;
		}

		for (int i=0; i < nT; i++){
			t(i + 1) = i * dt;
			RK4(i, t(i));
		}
	}

	void RK4(int i, double t){
		double h = dt / 2;
		arma::cube u, k1, k2, k3, k4;
		u = arma::cube(3, N, 2);  // 3dim x N paricles x 2 phases (postion and velocity)
		u.slice(0) = r.row(i);  // postion at time i
		u.slice(1) = v;  // velocity

		k1 = advance(t, u);
		k2 = advance(t + h, u + h * k1);
		k3 = advance(t + h, u + h * k2);
		k4 = advance(t + 2 * h, u + h * k3);
		u += h * (k1 + 2 * k2 + 2 * k3 + k4) / 3;

		r.row(i + 1) = u.slice(0);
		v = u.slice(1);
	}

	arma::cube advance(double time, arma::cube u){
		arma::cube du(size(u));
		arma::mat F = sum_particle_forces(u.slice(0));  // calculate ppi from poistions
		double V = get_Efield_at_time(time);
		du.slice(0) = u.slice(1);  // set change in pos to vel
		du.slice(1).row(0) = w0 * u.slice(1).row(1) + wzsq * u.slice(0).row(0) * V0 / 2 + F.col(0);
		du.slice(1).row(1) = -w0 * u.slice(1).row(0) + wzsq * u.slice(0).row(1) * V0 / 2 + F.col(1);
		du.slice(1).row(2) = -wzsq * u.slice(0).row(2) * V0 + F.col(2);

		return du;
	}


};

class Solver{
  public:
	int nT;
	double T, dt;
	arma::cube U;
	arma::vec t;

	Solver(double Time, double timestep, arma::mat U0, bool Rk4=true){
		dt = timestep;
		nT = (int)(Time / dt);
		t = arma::vec(nT, arma::fill::zeros);
		U = arma::cube(nT, U0.n_rows, U0.n_cols);
		U.zeros();
		U.row(0) = U0;
	}

	void solve(void(*func)(arma::mat)){
		for (int i=0; i < nT; i++){
			t(i + 1) = i * dt;
			RK4((*func), t(i + 1));
		}
	}
	void RK4(void(*func)(arma::mat), double t){

	}

};

double f(double t){
	return sin(t);
}

int main() {
	// Particle p1 = Particle(arma::vec(3).fill(1), arma::vec(3).fill(2), 0.1, 3);
	// Particle p2 = Particle(arma::vec(3).fill(3), arma::vec(3).fill(2), 0.1, 3);
	// double b = 1.;
	// double v = 3.;
	// double d = 4.;
	// vector<Particle> p;
	// p.push_back(p1);
	// p.push_back(p2);

    // PenningTrap P = PenningTrap(b, (*f), d, true);
	// P.insert_particles(p);
	// P.simulate(2, 0.5);
	// // cout << P.get_Efield_at_time(0.2) << endl;
	// cout << P.sum_particles_forces() << endl;
	arma::mat U0(3, 2, arma::fill::randn);
	U0.print();
	arma::mat U(size(U0));
	U.print();
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
