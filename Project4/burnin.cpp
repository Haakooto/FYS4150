#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>

#include "utils.cpp"

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
	/*
	Study burnin time for L=20 for T=1 and T=2.4
	Uses the 3 different initializations, and saves e and m to file.
	*/

	int M, L = 20;
	string fname;
	if (argc != 3){
		std::cout << "Bad usage! This program takes two params";
		std::cout << "\nfilename, and number of monte carlo cycles\n";
		return 1;
	}else{
		fname = argv[1];
        M = atoi(argv[2]);
    }
	arma::mat data;

	// open outfile
	ofstream out;
	out.open(fname + ".csv");
	out << "er1 mr1 er24 mr24 el1 ml1 el24 ml24 eh1 mh1 eh24 mh24\n";

	// loop over initializations
	for (const char* start:{"random", "lowest", "highest"}){
		for (int T: {1., 2.4}){  // loop over temperatures
			// run cycles
			arma::mat run = mc_run_culm(L, M, T, start, 0);
			// get data
			arma::vec e = run.row(0).t();
			arma::vec m = run.row(2).t();
			data = arma::join_rows(data, e);
			data = arma::join_rows(data, m);
		}
	}
	data.save(out, arma::raw_ascii);
	return 0;
}