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
	Dette er kun en kopi av burning.cpp
	For å gjøre den til tempting, endre input parametre
	og loop over temperatur i stedenfor initializeringen
	Legg til pragma og kall multi_mc, ikke mc_run_cuml
	Ellers skal alt være blomster og puter, mer eller mindre
	*/

	int M, T, L = 20;
	string fname;
	if (argc != 4){
		cout << "Bad usage! This program takes three params";
		cout << "\nfilename, temperature, and number of monte carlo cycles\n";
		return 1;
	} else {
		fname = argv[1];
		T = atoi(argv[2]);
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

		// get data
		data = join_rows(data, run.col(0));  // energy
		data = join_rows(data, run.col(1));  // energy, cumulative average
		data = join_rows(data, run.col(2));  // energy
		data = join_rows(data, run.col(3));  // energy, cumulative average
	}

	data.save(out, csv_ascii);
	return 0;
}