#include "armadillo"

arma::mat make_A(int, double);
void write_to_file(int, int, arma::mat&, arma::vec&, arma::mat&, arma::vec&);

arma::vec analytic_eigvec(arma::vec&, arma::mat&, int, double, double);
double analytic_eigval(int, int, double, double);
void analytic_solutions(arma::vec&, arma::mat&, int, double, double);

double max_diff(arma::mat, arma::mat);
double max_diff(arma::vec, arma::vec);
void sort_mat_by_vec(arma::mat&, arma::vec&);

