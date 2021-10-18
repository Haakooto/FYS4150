#include "PenningTrap.hpp"
extern const double ke;


PenningTrap::PenningTrap(double Bfield, double Efield, double length, bool particle_particle){
	B0 = Bfield;
	V0 = Efield;
	d = length;
	time_dep_V = false;
	ppi = particle_particle;
}

PenningTrap::PenningTrap(double Bfield, double(*tEfield)(double), double length, bool particle_particle){
	B0 = Bfield;
	tV0 = tEfield;
	d = length;
	time_dep_V = true;
	ppi = particle_particle;
}

void PenningTrap::insert_particles(vector<Particle> P){
	for (Particle p: P){
		particles.push_back(p);
		N++;
	}
}

void PenningTrap::insert_particles(Particle p){
	particles.push_back(p);
	N++;
}

void PenningTrap::insert_particles(int n, double m, int q){ // inserting n identical particles uniformly in sphere
	// https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
	for (int i=0; i<n; i++){
		arma::vec r(3, arma::fill::randn);
		arma::vec u(1, arma::fill::randu);
		r *= d * pow(u(0), 3) / arma::norm(r);
		particles.push_back(Particle(r, arma::zeros(3), m, q));
		N++;
	}
}

double PenningTrap::get_Efield_at_time(double t){
	if (time_dep_V) {return tV0(t);}
	else {return V0;}
}

arma::mat PenningTrap::sum_particle_forces(arma::mat ri){
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

void PenningTrap::analytic(double T, double timestep, double x0, double z0, double y_v0){   //analytic solution for one particle with a given starting position
	dt = timestep;
	nT = (int)(T / dt) + 1;
	// t = arma::vec(nT, arma::fill::zeros);
	r_a = arma::mat(nT, 3);  // nT timesteps x 3dim, only 1 particle
	r_a.row(0) = arma::vec({x0, 0, z0}).t();   //initial position

	double w_0 = particles[0].q * B0 / particles[0].m;
	double w_z = 2*particles[0].q * V0/(particles[0].m * pow(d, 2));

	complex <double> w_plus ((w_0+pow(pow(w_0, 2) - 2 * w_z, 0.5))/2, 0);
	complex <double> w_min (((w_0)-pow(pow(w_0, 2) - 2 * w_z, 0.5))/2, 0);

	complex <double> A_plus = (y_v0 + w_min*x0)/(w_min - w_plus);
	complex <double> A_min = -(y_v0 + w_plus*x0)/(w_min - w_plus);

	double time = dt;

	complex <double> complex_i(0, 1);

	for (int i=1; i < nT; i++){

		complex <double> f = A_plus * exp(- complex_i * w_plus * (complex<double>) time) + A_min * exp(-complex_i * w_min * (complex<double>)time);
		double x = real(f);
		double y = imag(f);
		double z = z0 * cos(pow(real(w_z), 0.5)*time);

		r_a.row(i) = arma::vec({x,y,z}).t();
		time += dt;
	}

}

void PenningTrap::simulate(double T, double timestep, string method){
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
        if (method == "RK4"){
            RK4(i, t(i));
        }
        else if (method == "Euler"){
            Euler(i, t(i));
        }

		t(i + 1) = (i + 1) * dt;

        for (int p=0; p < N; p++){
            if (arma::norm(r.slice(p).row(i + 1)) > d){
                Q(p) = 0;
            }
        }
	}
}


void PenningTrap::Euler(int i, double t){
    double h = dt;
    arma::cube u;
    u = arma::cube(3, N, 2);  // 3dim x N particles x 2 phases (position and velocity)

	arma::mat r_(r.row(i));  // position at time i
	if (N == 1){ u.slice(0) = r_.t();}  // wacky bugfix.  For 1 particle, need to transpose r.row
	else { u.slice(0) = r_;}
	u.slice(1) = v;  // velocity

    u += h * advance(t, u);

    r.row(i + 1) = u.slice(0);
	v = u.slice(1);
}


void PenningTrap::RK4(int i, double t){
	double h = dt / 2;   //divide by 2 to get rid of factor of 1/2 in later expression
	arma::cube u, k1, k2, k3, k4;
	u = arma::cube(3, N, 2);  // 3dim x N particles x 2 phases (position and velocity)

	arma::mat r_(r.row(i));  // position at time i
	// wacky bugfix.  For 1 particle, need to transpose r.row
	if (N == 1){ u.slice(0) = r_.t();}
	else { u.slice(0) = r_;}
	u.slice(1) = v;  // velocity

	k1 = advance(t, u);
	k2 = advance(t + h, u + h * k1);
	k3 = advance(t + h, u + h * k2);
	k4 = advance(t + 2 * h, u + h * k3);
	u += h * (k1 + 2 * k2 + 2 * k3 + k4) / 3;

	r.row(i + 1) = u.slice(0);
	v = u.slice(1);
}


arma::cube PenningTrap::advance(double time, arma::cube u){
	arma::cube du(size(u));  // change in u
	arma::mat F = sum_particle_forces(u.slice(0));  // calculate ppi from poistions
	du.slice(0) = u.slice(1);  // set change in pos to vel

	Efield(u.slice(0), time);
	Bfield(u.slice(1));

	//update velocities for x, y, z for all particles
	du.slice(1).row(0) = u.slice(0).row(0) - u.slice(1).row(1) + F.row(0);
	du.slice(1).row(1) = u.slice(0).row(1) - u.slice(1).row(0) + F.row(1);
	du.slice(1).row(2) = u.slice(0).row(2) + F.row(2);

	return du;  // 3dim x N particles x 2 phases (position and velocity)
}

void PenningTrap::Efield(arma::mat &pos, double t){
	double V = get_Efield_at_time(t);
	pos.row(0) %= Q / M * pow(d, -2) * V;
	pos.row(1) %= Q / M * pow(d, -2) * V;
	pos.row(2) %= -2 * Q / M * pow(d, -2) * V;
}

void PenningTrap::Bfield(arma::mat &vel){
	vel.row(0) %= Q * B0 / M;
	vel.row(1) %= -Q * B0 / M;
}

arma::cube PenningTrap::get_history(){
	return r;
}

arma::mat PenningTrap::get_asol(){
	return r_a;
}

arma::vec PenningTrap::get_time(){
	return t;
}
int PenningTrap::escaped(){
	int out = 0;
	for (int i=0; i < N; i++){
		if (particles[i].q != Q(i)){
			out++;
		}
	}
	return out;
}
