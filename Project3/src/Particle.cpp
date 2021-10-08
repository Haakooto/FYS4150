#include "Particle.hpp"

Particle::Particle(double charge, double mass, arma::vec position, arma::vec velocity){
    q = charge;
    m = mass;
    r = position;
    v = velocity;
}
