#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>
#include <omp.h>
#include <chrono>

#include "ising_model.cpp"

using namespace std;
using namespace arma;
using namespace std::chrono;

int main(int argc, char* argv[]) {
	/*
    Measures the time spent running the algorithm. 

    Arguments:
        M: int
            No. of Monte Carlo cycles
        R: int
            Number of concurrent runs
        L: int
            Lattice size
        T: float
            Temperature
        para: bool
            wether to parallelize loop over R
	*/
    bool para;
    int M, R, L;
    float T;

	string fname;
	if (argc != 6){
		cout << "Bad usage! This program takes five params";
		cout << "\n MC cycles, Runs, L, temperature, parabool \n";
		return 1;
	} else {
        M = atoi(argv[1]);
        R = atoi(argv[2]);
        L = atoi(argv[3]);
        T = atof(argv[4]);
        para = string(argv[5]) == "para";
    }

    string method = "random";
    int burnin = 100;
    arma::vec run(8, arma::fill::zeros);

    steady_clock::time_point t1 = steady_clock::now();
    multi_mc(L, M, R, T, run, method, burnin, para);
    steady_clock::time_point t2 = steady_clock::now();

    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    cout << time_span.count() << endl;

	return 0;
}
