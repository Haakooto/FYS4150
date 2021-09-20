#include "utils.hpp"

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
