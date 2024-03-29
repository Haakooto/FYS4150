#include <armadillo>
#include <functional>

using namespace std;

#ifndef __Particle_hpp__
#define __Particle_hpp__

class Particle
{
public:
    int q;
    double m;
    arma::vec r, v;

    Particle(arma::vec position, arma::vec velocity, double mass, int charge);
    void print();
};

#endif


#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

class PenningTrap
{
public:
    int N=0;  // number of particles in trap
    PenningTrap(double, double, double, bool=true);
    PenningTrap(double, function<double(double)>, double, bool=false);

    // Overloaded insertion methods
    void insert_particles(vector<Particle> P);
    void insert_particles(Particle);
    void insert_particles(int, double, int);

    void set_tEfield(function<double(double)>);
    void simulate(double, double, string="RK4");
    arma::vec analytic_analysis(double, double, string);

    arma::cube get_history();  // returns r and v
    arma::mat get_asol();  // returns r_a
    arma::vec get_time();  // retruns t
    int escaped();  // counts particles left in trap
    bool ppi;

private:
    int nT;  // number of time steps
    double dt;  // length of time step
    double B0, V0, d;  // properties of trap
    double r_cutoff, cut=0.3;  // cutoff for coulomb force, set to 0.3d
    function<double(double)> tV0; // time-dep E-field (lambda func)
    vector<Particle> particles;  // Particle container
    bool time_dep_V;
    arma::cube R;  // postion and velocity of all particles at all times
    arma::mat r_a; // analytic solutions.
    arma::vec t;  // time vector
    arma::rowvec Q, M;  // properties of particles

    double get_Efield_at_time(double);
    arma::mat sum_particle_forces(arma::mat);
    void RK4(arma::cube &, double);
    void Euler(arma::cube &, double);
    arma::cube advance(double, arma::cube);
    void Efield(arma::mat &, double );
    void Bfield(arma::mat &);
    void analytic_sols(double, double, double);

};

#endif
