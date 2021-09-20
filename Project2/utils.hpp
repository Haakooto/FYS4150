#include "armadillo"

using namespace arma;

mat make_A(int, double);

vec analytic_eigvec(vec&, mat&, int, double, double);
double analytic_eigval(int, int, double, double);
void analytic_solutions(vec&, mat&, int, double, double);

double max_diff(mat, mat);
double max_diff(vec, vec);
void sort_mat_by_vec(mat&, vec&);
