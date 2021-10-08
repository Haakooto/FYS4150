#include <iostream>
#include <armadillo>

#include "Particle.hpp"

int main(){
    Particle p1(3.0, 2.0, arma::vec(3).fill(2.), arma::vec(1.).fill(4.));
    std::cout << p1.q << p1.m << p1.r << p1.v << std::endl;

    return 0;
}