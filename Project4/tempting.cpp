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
    Run MC-cylces over many temperatures

	Input: filename, M, R, L, Tmin, Tmax, Ts
	Output: csv-file with measured quantities for each temperature
    */
    int M, L, Ts;   //Ts is the number of temperatures to iterate over
    double Tmin, Tmax;

	string fname;
	if (argc != 7){
		cout << "Bad usage! This program takes six params";
		cout << "\n filename, MC cycles, L, Tmin, Tmax, # temperatures \n";
		return 1;
	} else {
        fname = argv[1];
        M = atoi(argv[2]);
        L = atoi(argv[3]);
        Tmin = atof(argv[4]);
        Tmax = atof(argv[5]);
        Ts = atoi(argv[6]);
    }
	mat data(Ts, 9);

    double inc = (Tmax - Tmin)/(Ts - 1);
    string method = "random";
    int burnin = 500;

    #pragma omp parallel for
    for (int i = 0; i < Ts; i++)
    {
        double T = Tmin + inc * i;
        arma::vec run(8, arma::fill::zeros);
        single_mc(L, M, T, run, method, burnin);
        data(i, 0) = T;
        data.submat(i, 1, i, 8) = run.t();
    }

	// open outfile
	ofstream out;
	out.open("data/" + fname + ".csv");
    out << "L = " << L << " MC cycles = " << M << " Repetitions = " << R << endl;
	out << "T,e_avg,m_avg,Cv,chi,e_err,m_err,Cv_err,chi_err" << endl;
	data.save(out, csv_ascii);
	return 0;
}
