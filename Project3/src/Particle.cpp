#include "Particle.hpp"

Particle::Particle(arma::vec position, arma::vec velocity, double mass, int charge){
		r = position;
		v = velocity;
		m = mass;
		q = charge;
}
Particle::print(){
		cout << "Particle at ( ";
		for (int i=0; i<3; i++){cout << r(i) << " ";}
		cout << ") with velocity ( ";
		for (int i=0; i<3; i++){cout << v(i) << " ";}
		cout << ") with mass " << m << " and charge " << q << endl;
}
