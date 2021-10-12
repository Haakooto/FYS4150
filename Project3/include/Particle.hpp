#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle
{
public:
    double q, m;
    arma::vec r, v;

    Particle(double charge, double mass, arma::vec position, arma::vec velocity);
    arma::vec force_from_other(Particle other);
};

#endif