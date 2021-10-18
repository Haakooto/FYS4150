#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <vector>
#include <complex>

class PenningTrap
/*
*/
{
public:
    int N=0;  // number of particles in trap
    PenningTrap(double Bfield, double Efield, double length, bool particle_particle=false );
    PenningTrap(double Bfield, double(*)(double), double length, bool particle_particle=false );
    // Overloaded insertion methods
    void insert_particles(vector<Particle>);
    void insert_particles(Particle);
    void insert_particles(int, double, int);
    void simulate(double, double);
    void analytic(double, double, double, double, double);
    arma::cube get_history();  // returns r
    arma::vec get_time();  // retruns t
    int escaped();  // counts particles left in trap

private:
    int nT;  // number of time steps
    double dt;  // length of time step
    double B0, V0, d;  // properties of trap
    double (*tV0)(double);  // time-dep E-field
    vector<Particle> particles;  // Particle container
    bool ppi, time_dep_V;
    arma::cube r;  // postion of all particles at all times
    arma::mat v;  // velocity of all particles, only one time step
    arma::vec t;  // time vector
    arma::rowvec Q, M;  // properties of particles

    double get_Efield_at_time(double);
    arma::mat sum_particle_forces(arma::mat);
    void RK4(int, double);
    void Euler(int, double);
    arma::cube advance(double, arma::cube);
    void Efield(arma::mat &, double );
    void Bfield(arma::mat &);

};

#endif