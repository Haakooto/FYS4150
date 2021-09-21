#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>
#include "utils.hpp"

using namespace arma;
using namespace std;

const double tolerance = 1e-4;

double max_offdiag_symmetric(const mat &A, int &k, int &l){
    int size = A.n_rows;
    double max = 0;
    //   i is rows and j is columns (maybe)
    for (int i=0; i < size; i++){
        for (int j=i + 1; j < size; j++){
            if (abs(A(i,j)) > max){
                max = abs(A(i,j));
                k = i;
                l = j;
            }
        }
    }
    return max;
}

void Rotation(mat &A, mat &R, int &k, int &l, int &size){
    double tau = (A(l,l) - A(k,k)) / (2 * A(k,l));
    double t;
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

bool Jacobi(mat &A, mat &R, const double tol, const int maxiter, int &iters){
    int size = A.n_rows;
    int k;
    int l;
    double max = max_offdiag_symmetric(A, k, l);
    while (max > tol){
        Rotation(A, R, k, l, size);
        max = max_offdiag_symmetric(A, k, l);
        iters++;
        if (iters == maxiter){
            return false;
        }
    }
    return true;
}

void test_max_offdiag(){
    // Problem 4b
    mat A(4, 4, fill::eye);
    A(0, 3) = A(3,0) = 0.5;
    A(1, 2) = A(2,1) = -0.7;
    int k;
    int l;
    double max = max_offdiag_symmetric(A, k, l);
    assert (max == 0.7);
    assert (k == 1);
    assert (l == 2);
}

void test_analyticity(){
    // Problem 5b
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

void run_Jacobi(int N, const int maxiter){
	double hsq = (double)pow(N + 1, 2);  // h^-2
    int iters = 0;
	
	mat A = make_A(N, hsq);
    mat eig_vecs = mat(N, N, fill::eye);

	if (Jacobi(A, eig_vecs, tolerance, maxiter, iters)){
        cout << "Algorithm converged after " << iters << " iterations.\n";
        cout << "Writing results to file.\n";
        
        vec eig_vals = A.diag();
        sort_mat_by_vec(eig_vecs, eig_vals);


	    int vecs_to_write = 3;
        // calculate analytic solutions
	    vec a_vals(vecs_to_write);
	    mat a_vecs(N, vecs_to_write);
	    analytic_solutions(a_vals, a_vecs, N, 2 * hsq, -hsq);

        write_to_file(N, vecs_to_write, iters, eig_vecs, eig_vals, a_vecs, a_vals);     

    } else {
        cout << "Algorithm did not converge. Terminated after " << maxiter << " iterations.\n";
    }
}

int main() {
    tests();  // Run all test functions

	int N = 101;
    int maxiter = 100000; // Set to huge value. Should be estimated
    for (int i = 5; i < N ; i++){
        run_Jacobi(i, maxiter);
    }

	return 0;
}