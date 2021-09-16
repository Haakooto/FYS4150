#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;


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

int main() {
	int N = 6;
	int n = N + 1;
	double hsq = (double)pow(n, 2);

	mat A = make_A(N, hsq);
	vec D;
	mat S;

	eig_sym(D, S, A);

	int number_of_analytic_eig_to_check = 3;
	vec A_vals(number_of_analytic_eig_to_check);
	mat A_vecs(N, number_of_analytic_eig_to_check);

	analytic_solutions(A_vals, A_vecs, N, A(0, 0), A(0, 1));

	// A.print();
	D.print();
	S.print();

	A_vals.print();
	A_vecs.print();

	return 0;
}