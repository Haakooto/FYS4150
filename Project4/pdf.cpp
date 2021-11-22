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
    /*
    Generates the data for the estimated probability distributions of the energies for temperatures T=1 and T=2.4
    by calling on the multi_prob() function in ising_model. Writes the results to files.

    Arguments:
        M: int
            No. of Monte Carlo cycles
        R: int
            No. of repetitions for each cycle
    */

    int M, R, L=20;
    int burnin = 100;
	if (argc != 3){
		std::cout << "Bad usage! This program takes four params";
		std::cout << "\n number of monte carlo cycles, and number of cocurrent mcc\n";
		return 1;
	}else{
        M = atoi(argv[1]);
        R = atoi(argv[2]);
    }

    ofstream out1, out2;
    out1.open("data/pdf_T1.csv");
    out1 << "e_avg,prob" << endl;

    out2.open("data/pdf_T2.4.csv");
    out2 << "e_avg,prob" << endl;

    mat Prob = multi_prob(L, M, R, 1, burnin);
    Prob.save(out1, csv_ascii);

    Prob = multi_prob(L, M, R, 2.4, burnin);
    Prob.save(out2, csv_ascii);

    return 0;
}
