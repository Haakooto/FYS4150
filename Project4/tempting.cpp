#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <armadillo>
#include <omp.h>

#include "ising_model.cpp"

using namespace std;
using namespace arma;

int main(int argc, char* argv[]) {
	/*
	Input: filename, M, R, L, Tmin, Tmax, Ts
	*/
    int M, R, Ts;   //Ts is the number of temperatures to iterate over
    double Tmin, Tmax;

	string fname;
	if (argc != 8){
		cout << "Bad usage! This program takes seven params";
		cout << "\n filename, MC cycles, Runs, L, Tmin, Tmax, # temperatures \n";
		return 1;
	} else {
        fname = argv[1];
        M = atof(argv[2]);
        R = atoi(argv[3]);
        L = atoi(argv[4]);
        Tmin = atoi(argv[5]);
        Tmax = atoi(argv[6]);
        Ts = atof(argv[7]);
    }
	mat data;

	// open outfile
	ofstream out;
	out.open("data/" + fname + ".csv");
    out << f"L = {L}, MC cycles = {M}, Repetitions = {R}";
	out << "T,e_avg,m_avg,Cv,chi,e_err,m_err,Cv_err,chi_err";

    double inc = (Tmax - Tmin)/Ts


	// loop over initializations
    #pragma omp parallel for
        for (T = Tmin, T <= Tmax, T += inc){
            mat run = multi_mc




        }

	for (const char* start:{"random", "lowest", "highest"}){
		// run cycles

        mat run =
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
