#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>
#include <random>

#include "ising_model.cpp"

using namespace std;
using namespace arma;

int main(int argc, char* argv[]){
    
    string fname1, fname2;
    int M;
	if (argc != 4){
		std::cout << "Bad usage! This program takes three params";
		std::cout << "\nfilename1, filename2, and number of monte carlo cycles\n";
		return 1;
	}else{
		fname1 = argv[1];
        fname2 = argv[2];
        M = atoi(argv[3]);
    }
    ofstream out1, out2;
    out1.open("data/" + fname1 + ".csv");
    out2.open("data/" + fname2 + ".csv");
    
    mat Lattice = make_sys(20, "random");
    mat Prob = mc_run(Lattice, 1., M);
    Prob.save(out1, raw_ascii);
    
    Prob = mc_e_prob(Lattice, 2.4, M);
    Prob.save(out2, raw_ascii);

    return 0;
}