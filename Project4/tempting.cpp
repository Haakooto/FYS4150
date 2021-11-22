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
    Generates the data needed to find the critical temperature by calculating all the needed variables as
    function of temperature for the given range by calling the single_mc() function from ising_model.
    The results are written to a file.

    Arguments:
        fname: string
            Desired filename
        M: int
            No. of Monte Carlo cycles
        R: int
            No. of repetitions for each cycle
        L: int
            Lattice size
        Tmin: double
            Initial temperature
        Tmax: double
            Final temperature
        Ts: int
            No. of temperatures to iterate over
    */

    int M, R, L, Ts;
    double Tmin, Tmax;

	string fname;
	if (argc != 8){
		cout << "Bad usage! This program takes seven params";
		cout << "\n filename, MC cycles, Runs, L, Tmin, Tmax, # temperatures \n";
		return 1;
	} else {
        fname = argv[1];
        M = atoi(argv[2]);
        R = atoi(argv[3]);
        L = atoi(argv[4]);
        Tmin = atof(argv[5]);
        Tmax = atof(argv[6]);
        Ts = atoi(argv[7]);
    }
	mat data(Ts, 9);

	// open outfile
	ofstream out;
	out.open("data/" + fname + ".csv");
    out << "L = " << L << " MC cycles = " << M << " Repetitions = " << R << endl;
	out << "T,e_avg,m_avg,Cv,chi,e_err,m_err,Cv_err,chi_err" << endl;

    double inc = (Tmax - Tmin)/(Ts - 1);
    int i = 0;
    string method = "random";
    int burnin = 1000;

    // omp_set_num_threads(4);
	// loop over initializations
    #pragma omp parallel for
        for (i = 0; i < Ts; i++)
        {
            double T = Tmin + inc * i;
            arma::vec run(8, arma::fill::zeros);
            single_mc(L, M, T, run, method, burnin);
            data(i, 0) = T;
            data.submat(i, 1, i, 8) = run.t();
        }

	data.save(out, csv_ascii);
	return 0;
}
