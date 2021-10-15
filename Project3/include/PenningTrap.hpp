#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <vector>

class PenningTrap
/*
*/
{
public:
    int N=0;
    PenningTrap(double Bfield, double Efield, double length, bool particle_particle=false );
    PenningTrap(double Bfield, double(*)(double), double length, bool particle_particle=false );
    void insert_particles(vector<Particle>);
    void insert_particles(Particle);
    void insert_particles(int, double, int);
    void simulate(double, double);
    void analytic(double, double, double, double, double);

private:
    int nT;
    double dt;
    double B0, V0, d;
    double (*tV0)(double);
    vector<Particle> particles;
    bool ppi, time_dep_V;
    arma::cube r;
    arma::mat v;
    arma::vec t;
    arma::rowvec Q, M;

    double get_Efield_at_time(double);
    arma::mat sum_particle_forces(arma::mat);
    void RK4(int, double);
    void Euler(int, double);
    arma::cube advance(double, arma::cube);
    void Efield(arma::mat &, double );
    void Bfield(arma::mat &);


};

#endif