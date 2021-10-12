#include "Particle.hpp"

Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity){
    q = charge;
    m = mass;
    r = position;
    v = velocity;
}

Particle::force_from_other(Particle other){
    rirj = r - other.r
    norm = arma::norm(rirj)
    return other.q * rirj * pow(norm, -3)
}