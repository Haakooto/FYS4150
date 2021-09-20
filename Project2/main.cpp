#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>
#include "armadillo"
#include "utils.hpp"

using namespace std;
using namespace arma;

const double tolerance = 1e-14;

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
    return abs(max);
}

void Rotation(mat &A, mat &R, int &k, int &l, int &size){
    double tau = (A(l,l) - A(k,k)) / (2 * A(k,l));
    double t;
    // Hvis Håkon Olav kjenner en enkel sign function kan vi endre her
    if (tau > 0){
        t = 1 / (tau + sqrt(1 + pow(tau, 2)));
    } else {
        t = - 1 / (- tau + sqrt(1 + pow(tau, 2)));
    }

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
            A(l,i) = A(i,l);
        }
        r_ik = R(i,k);
        R(i,k) = R(i,k) * c - R(i,l) * s;
        R(i,l) = R(i,l) * c + r_ik * s;
    }
}

mat Jacobi(mat &A, double tol){
    int size = A.n_rows;
    int k;
    int l;
    int counter = 1;
    double max = max_offdiag_symmetric(A, k, l);
    mat R = mat(size, size, fill::eye);
    while (max > tol){
        Rotation(A, R, k, l, size);
        max = max_offdiag_symmetric(A, k, l);
    }
    return R;
}

void test_max_offdiag(){
    mat A(4, 4, fill::eye);
    A(0, 3) = A(3,0) = 0.5;
    A(1, 2) = A(2,1) = - 0.7;
    int k;
    int l;
    double max = max_offdiag_symmetric(A, k, l);
    assert (max == 0.7);
    assert (k == 1);
    assert (l == 2);
}

void test_analyticity(){
    int N = 6;
    double hsq = (double)pow(N + 1, 2);
    double tol = 1e-10;

    // using armadillo to solve problem
    mat A = make_A(N, hsq);
    vec D; mat S;
    eig_sym(D, S, A);

    // finding analytic solutions
    vec avals(N);
    mat avecs(N, N);
    analytic_solutions(avals, avecs, N, 2 * hsq, -hsq);

    // comparing the two
    sort_mat_by_vec(S, D);
    sort_mat_by_vec(avecs, avals);
    assert (max_diff(S, avecs) < tol);
    assert (max_diff(D, avals) < tol);
}

void tests(){
    test_analyticity(); // problem 3 in project discription
    test_max_offdiag();

    cout << "All tests passed" << endl;
}

int main() {
    tests();  // Run all test functions

	int N = 6;
	int n = N + 1;
	int k = 0;
	int l = 0;
	double hsq = (double)pow(n, 2);
	

	mat A = make_A(N, hsq);
	vec D;
	mat S;

	// eig_sym(D, S, A);

	int number_of_analytic_eig_to_check = N;
	vec a_vals(number_of_analytic_eig_to_check);
	mat a_vecs(N, number_of_analytic_eig_to_check);


	analytic_solutions(a_vals, a_vecs, N, A(0, 0), A(0, 1));

	mat R = Jacobi(A, tolerance);
    vec v = A.diag();
    sort_mat_by_vec(R, v);

    
	cout << " " << endl;
	R.print();
	cout << " " << endl;
    // a_vecs.print();
    v.print();
    // cout << max_diff(R, a_vecs) << endl;
	
    // // A.print();
	// // D.print();
	// // S.print();

	// a_vals.print();
	// a_vecs.print();

	// max_test();

	return 0;
}