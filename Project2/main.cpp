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

double max_offdiag_symmetric(const mat &A, int &k, int &l){
    int size = A.n_rows;
    double max = 0;
    //   i is rows and j is columns (maybe)
    for (int i=0; i < size; i++){
        for (int j=i + 1; j < size; j++){
            if (abs(A(i,j)) > abs(max)){
                max = A(i,j);
                k = i;
                l = j;
                }
            }
        }
    return max;}

void max_test(){
    mat A(4, 4, fill::eye);
    A(0, 3) = A(3,0) = 0.5;
    A(1, 2) = A(2,1) = - 0.7;
    int k;
    int l;
    double max = max_offdiag_symmetric(A, k, l);
    cout << max << endl << k << endl << l << endl;
    }

void Rotation(mat &A, mat &R, int &k,int &l, int &size){
    double tau = (A(l,l) - A(k,k)) / (2 * A(k,l));
    double t;
    // Hvis Håkon Olav kjenner en enkel sign function kan vi endre her
    if (tau > 0){
        t = 1 / (tau + sqrt(1 + pow(tau, 2)));}
    else if (tau < 0){
        t = - 1 / (- tau + sqrt(1 + pow(tau, 2)));}
    else{t = 1;};
    double c = 1 / sqrt(1 + pow(t, 2));
    double s = c * t;

    double a_kk = A(k,k);
    A(k,k) = A(k,k) * pow(c, 2) - 2 * A(k,l) * c * s + A(l,l) * pow(s, 2);
    A(l,l) = A(l,l) * pow(c, 2) + 2 * A(k,l) * c * s + a_kk * pow(s, 2);
    A(k,l) = A(l,k) = 0;



    double r_ik;

    for (int i=0; i < size; i++){
        if ((i != k) && (i != l)){
            A(i,k) = A(i,k) * c - A(i,l) * s;
            A(i,l) = A(i,l) * c + A(k,i) * s;
            A(k,i) = A(i,k);
            A(l,i) = A(i,l);};

        r_ik = R(i,k);
        R(i,k) = R(i,k) * c - R(i,l) * s;
        R(i,l) = R(i,l) * c + r_ik * s;
        };
    }

mat Jacobi(mat &A, double tol){
    int size = A.n_rows;
    int k;
    int l;
    int counter = 1;
    double max = max_offdiag_symmetric(A, k, l);
    mat R = mat(size, size, fill::eye);
    while (abs(max) > tol){
        Rotation(A, R, k, l, size);
        max = max_offdiag_symmetric(A, k, l);
        };
    return R;
    }

int main() {
	int N = 6;
	int n = N + 1;
	int k = 0;
	int l = 0;
	double hsq = (double)pow(n, 2);
	double tolerance = 1e-8;

	mat A = make_A(N, hsq);
	vec D;
	mat S;

	// eig_sym(D, S, A);

	int number_of_analytic_eig_to_check = 3;
	vec A_vals(number_of_analytic_eig_to_check);
	mat A_vecs(N, number_of_analytic_eig_to_check);



	analytic_solutions(A_vals, A_vecs, N, A(0, 0), A(0, 1));

	mat R = Jacobi(A, tolerance);

	R.print();
	cout << " " << endl;
	A.print();
	// A.print();
	// D.print();
	// S.print();

	A_vals.print();
	A_vecs.print();

	// max_test();

	return 0;
}