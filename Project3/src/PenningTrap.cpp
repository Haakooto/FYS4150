#include "PenningTrap.hpp"
extern const double ke;


PenningTrap::PenningTrap(double Bfield, double Efield, double length, bool particle_particle){
	B0 = Bfield;
	V0 = Efield;
	d = length;
	time_dep_V = false;
	ppi = particle_particle;
}

PenningTrap::PenningTrap(double Bfield, function<double(double)> Efield, double length, bool particle_particle){
	B0 = Bfield;
	tV0 = Efield;
	d = length;
	time_dep_V = true;
	ppi = particle_particle;
}

void PenningTrap::insert_particles(vector<Particle> P){
	// insert many particles
	for (Particle p: P){
		particles.push_back(p);
		N++;
	}
}

void PenningTrap::insert_particles(Particle p){
	// Insert single particle
	particles.push_back(p);
	N++;
}

void PenningTrap::insert_particles(int n, double m, int q){
	// Insert n identical particles with random position and velocity
	for (int i=0; i<n; i++){
        arma::vec r = arma::vec(3).randn() * 0.1 * d;  // random initial position
        arma::vec v = arma::vec(3).randn() * 0.1 * d;  // random initial velocity

		particles.push_back(Particle(r, v, m, q));
		N++;
	}
}

void PenningTrap::set_tEfield(function<double(double)> f){
	// change function for time-dep Efield
	tV0 = f;
}

double PenningTrap::get_Efield_at_time(double t){
	// evaluate Efield at time t. May be constant
	if (time_dep_V) {return tV0(t);}
	else {return V0;}
}

void PenningTrap::simulate(double T, double timestep, string method){
	dt = timestep;
	nT = (int)(T / dt) + 1;
	t = arma::vec(nT, arma::fill::zeros);
	R = arma::cube(6, N, nT);  // (3dim x 2 phases) : N particles : nT timesteps

	Q = arma::rowvec(N);  // charges
	M = arma::rowvec(N);  // masses
	r_cutoff = cut * d;


	// initiate simulation, 
	for (int i=0; i < N; i++){
		Particle p = particles[i];
		R.slice(0).rows(0, 2).col(i) = p.r; // set positions
		R.slice(0).rows(3, 5).col(i) = p.v;  // set velocities

		Q(i) = p.q;  // set charges
		M(i) = p.m; // set masses
	}

	// start simulation
	arma::cube u(3, N, 2);  // 3dim : N particles : 2 phases
	for (int i=0; i < nT - 1; i++){
		u.slice(0) = R.slice(i).rows(0, 2);  // positions at time i
		u.slice(1) = R.slice(i).rows(3, 5);  // velocities at time i

		if (method == "Euler"){
			Euler(u, t(i));
		} else{
			RK4(u, t(i));
		}

		R.slice(i + 1).rows(0, 2) = u.slice(0);  // write new positions
		R.slice(i + 1).rows(3, 5) = u.slice(1);  // write new velocities
		t(i + 1) = (i + 1) * dt; 

		// Check if any particle is outside trap
        for (int p=0; p < N; p++){
            if (arma::norm(R.slice(i + 1).rows(0, 2).col(p)) > d){
                Q(p) = 0;  // effectively set E and B-field to 0 outside trap
            }
        }
	}
}

void PenningTrap::Euler(arma::cube &u, double t){
    u += dt * advance(t, u);
}

void PenningTrap::RK4(arma::cube &u, double t){
	double h = dt / 2;   //divide by 2 to get rid of factor of 1/2 in later expression
	arma::cube k1, k2, k3, k4;

	k1 = advance(t, u);
	k2 = advance(t + h, u + h * k1);
	k3 = advance(t + h, u + h * k2);
	k4 = advance(t + 2 * h, u + 2 * h * k3);
	u += h * (k1 + 2 * k2 + 2 * k3 + k4) / 3;
}


arma::cube PenningTrap::advance(double time, arma::cube u){
	arma::cube du(size(u));  // change in u
	arma::mat F = sum_particle_forces(u.slice(0));  // calculate ppi from poistions
	du.slice(0) = u.slice(1);  // set change in pos to vel

	Efield(u.slice(0), time);
	Bfield(u.slice(1));

	// change in velocities for x, y, z for all particles
	du.slice(1).row(0) = F.row(0) + u.slice(0).row(0) - u.slice(1).row(1);
	du.slice(1).row(1) = F.row(1) + u.slice(0).row(1) - u.slice(1).row(0);
	du.slice(1).row(2) = F.row(2) + u.slice(0).row(2);

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

arma::mat PenningTrap::sum_particle_forces(arma::mat ri){
	/*
	Calculates all particle-particle forces
	*/
	arma::mat Eforce(N, 3, arma::fill::zeros);  // N particles : 3dim
	if (!ppi) {return Eforce.t();}   // ignore particle-particle interactions

	else {
		for (int i=0; i < N; i++){
			for (int j=0; j < i; j++){  
				arma::vec diff_ri = ri.col(i)-ri.col(j);
				double norm = arma::norm(diff_ri);
				if (norm < r_cutoff){  // speed-up
					arma::vec F = ke * Q(i) * Q(j) * diff_ri * pow(norm, -3);
					Eforce.row(i) += F.t() / M(i);  // Fij = -Fji
					Eforce.row(j) -= F.t() / M(j);
				}
			}
		}
	} return Eforce.t();
}

arma::vec PenningTrap::analytic_analysis(double T, double timestep, string method){
	ppi = false;  // make sure analytic sols are valid
	simulate(T, timestep, method);  // Do numertic simulation
	analytic_sols(particles[0].r[0], particles[0].r[2], particles[0].v[1]);  // Do analytic ones

	arma::cube np = R.rows(0, 2);  // positions
	arma::mat r = np.col(0);  // first particle
	arma::vec rel_errs = arma::vec(nT);
	arma::vec abs_errs = arma::vec(nT);

	for (int i=1; i < nT; i++){
		double abs_err = arma::norm(r_a.col(i) - r.col(i));
		double rel_err = abs_err / arma::norm(r_a.col(i));
		rel_errs(i) = rel_err;
		abs_errs(i) = abs_err;
	}
	rel_errs(0) = arma::max(abs_errs);  // error in first step is 0 anyways, store max abs err there
	return rel_errs;
}

void PenningTrap::analytic_sols(double x0, double z0, double y_v0){   //analytic solution for one particle with a given starting position
	r_a = arma::mat(3, nT);  // 3dim : nT timesteps, only 1 particle
	r_a.col(0) = arma::vec({x0, 0, z0});   //initial position

	double w_0 = particles[0].q * B0 / particles[0].m;
	double w_z = 2 * particles[0].q * V0 / (particles[0].m * pow(d, 2));

	complex <double> w_p ((w_0 + pow( pow(w_0, 2) - 2 * w_z, 0.5 )) / 2, 0);
	complex <double> w_m ((w_0 - pow( pow(w_0, 2) - 2 * w_z, 0.5 )) / 2, 0);

	complex <double> A_plus = (y_v0 + w_m * x0) / (w_m - w_p);
	complex <double> A_min = -(y_v0 + w_p * x0) / (w_m - w_p);

	double time = dt;

	complex <double> complex_i(0, 1);

	for (int i=1; i < nT; i++){

		complex <double> f = A_plus * exp(-complex_i * w_p * (complex<double>)time) + A_min * exp(-complex_i * w_m * (complex<double>)time);
		double x = real(f);
		double y = imag(f);
		double z = z0 * cos(pow(w_z, 0.5) * time);

		r_a.col(i) = arma::vec({x, y, z});
		time += dt;
	}
}

arma::cube PenningTrap::get_history(){
	return R;
}

arma::mat PenningTrap::get_asol(){
	return r_a;
}

arma::vec PenningTrap::get_time(){
	return t;
}

int PenningTrap::escaped(){
	// Particle charge change if it left trap
	int out = 0;
	for (int i=0; i < N; i++){
		if (particles[i].q != Q(i)){
			out++;
		}
	}
	return out;
}
