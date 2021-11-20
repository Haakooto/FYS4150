#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>

#include "ising_model.cpp"

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
	/*
	Study burnin time for L=20 for T=1 and T=2.4
	Uses the 3 different initializations, and saves e and m to file.
	*/

	int M, T, L = 20;
	string fname;
	if (argc != 4){
		cout << "Bad usage! This program takes three params";
		cout << "\nfilename, temperature, and number of monte carlo cycles\n";
		return 1;
	} else {
		fname = argv[1];
		T = atof(argv[2]);
        M = atoi(argv[3]);
    }
	mat data;

	// open outfile
	ofstream out;
	out.open("data/" + fname + ".csv");
	out << "e_rnd,avg_e_rnd,m_rnd,avg_m_rnd,";
	out << "e_low,avg_e_low,m_low,avg_m_low,";
	out << "e_hig,avg_e_hig,m_hig,avg_m_hig\n";

	// loop over initializations
	for (const char* start:{"random", "lowest", "highest"}){
		// run cycles
		mat run = mc_run_cuml(L, M, T, start, 0).t();

		data = join_rows(data, run);
	}

	data.save(out, csv_ascii);
	return 0;
}
