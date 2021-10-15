#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

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