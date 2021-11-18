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
    int M, R, L=20;
    int burnin = 0;
	if (argc != 5){
		std::cout << "Bad usage! This program takes four params";
		std::cout << "\nfilename1, filename2, number of monte carlo cycles, and number of cocurrent mcc\n";
		return 1;
	}else{
		fname1 = argv[1];
        fname2 = argv[2];
        M = atoi(argv[3]);
        R = atoi(argv[4]);
    }

    ofstream out1, out2;
    out1.open("data/" + fname1 + ".csv");
    out2.open("data/" + fname2 + ".csv");

    mat Prob = multi_prob(L, M, R, 1, burnin);
    Prob.save(out1, csv_ascii);

    Prob = multi_prob(L, M, R, 2.4, burnin);
    Prob.save(out2, csv_ascii);

    return 0;
}