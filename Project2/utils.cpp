#include <fstream>
#include <iomanip>
#include "utils.hpp"

using namespace arma;
using namespace std;

mat make_A(int N, double hsq){
	mat A(N, N, fill::zeros);
	for (int i = 0; i < N; i++){
		A(i, i) = 2 * hsq;
		if (i != 0){
			A(i, i - 1) = -hsq;
			A(i - 1, i) = -hsq;
		}
	}
	return A;
}

double analytic_eigval(int i, int N, double d, double a){
	return d + 2 * a * cos(i * M_PI / (N + 1));
}

vec analytic_eigvec(int i, int N){
	vec tmp(N);
	for (int j = 1; j <= N; j++){
		tmp(j - 1) = sin(j * i * M_PI / (N + 1));
	}
	return normalise(tmp);
}

void analytic_solutions(vec &vals, mat &vecs, int N, double d, double a){
	int M = vals.size();
	for (int i=0; i < M; i++){
		vals(i) = analytic_eigval(i + 1, N, d, a);
		vecs.col(i) = analytic_eigvec(i + 1, N);
	}
}

double max_diff(mat A, mat B){
    return abs(A - B).max();
}

double max_diff(vec a, vec b){
    return abs(a - b).max();
}

void sort_mat_by_vec(mat &A, vec &b){
    int N = A.n_cols;
    uvec indices = sort_index(b);
    uvec cols = linspace<uvec>(0, N-1, N);
    A = A.submat(cols, indices);
    b = sort(b);

    for (int i=0; i<N; i++){
        if (A(0, i) < 0){
            A.col(i) *= -1;
        }
    }
}

void write_to_file(int N, int M, mat &vecs, vec &vals, mat &avec, vec &almb){
	double h = 1 / ((double)N + 1);
	ofstream out;
	out.open("results/data_" + to_string(N) + ".txt");
	// Make header of datafile. First col is x
	// Next M are analytic results, next M are computed results
	out << "x";
	string u = " analytic_";
	for (int i = 0; i < M * 2; i++){
		out << u + to_string(i % M + 1);
		if (i == M - 1) {
			u = " computed_";
		}
	}
	out << endl;
	out << fixed << setprecision(8);

	// write eigenvalues
	out << "0"; // value under x, to be discarded
	vec v = almb;
	for (int i = 0; i < M * 2; i++){
		out << " " << v(i % M);
		if (i == M - 1){
			v = vals;
		} 
	}
	out << endl;
	
	// write first boundary condition
	for (int i = 0; i < M * 2 + 1; i++){
		out << 0 << " ";
	}
	out << endl;

	for (int j=0; j < N; j++){
		out << (double)(j + 1) * h;  // xj
		mat v = avec;
		for (int i=0; i < M * 2; i++){
			out << " " << v(j, i % M);
			if (i == M - 1){
				v = vecs;
			}
		}
		out << endl;
	}

	// write other boundary condition
	out << (double)(N + 1) * h;
	for (int i = 0; i < M * 2; i++){
		out << " " << 0;
	}
	out << endl;
	out.close();
}