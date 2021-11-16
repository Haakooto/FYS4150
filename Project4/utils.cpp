#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <armadillo>
#include <random>



arma::mat make_sys(int L, std::string method="random"){
    /*
    Initialises a Ising model lattice.

    Arguments:
        L: int
            Size of the periodic lattice
        method: std::string
            Method used to initialise the lattice.
            Can be one of "random", "lowest" or "highest".
            Default is set as "random".
    Returns:
        Lattice: arma::mat
            Object describing the lattice of spins.
    */

    arma::mat Lattice(L, L);
    if (method == "random"){
        /*
        Her må settes inn en annen rng, denne defaulten gir alltid de samme tallene.
        HO, jeg stoler på at du finner ut av det
        */
        arma_rng::set_seed_random();
        Lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0, 1));
        Lattice *= 2;
        Lattice -= 1;
    }
    else if (method == "lowest")
    {
        Lattice.ones();
    }
    else if (method == "highest"){
        Lattice.ones();
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                Lattice(i, j) *= pow(-1, i + j);
            };
        };
    }
    else{
        throw std::invalid_argument("invalid method");
    };

    return Lattice;
}

double calc_E(const arma::mat& Lattice){
    /*
    Calculates the initial energy of a system.

    Arguments:
        Lattice: arma::mat
            Object describing the lattice of spins.
    Returns:
        E: double
            The calculated energy of the system.
    */

    int R = Lattice.n_rows;
    int C = Lattice.n_cols;
    double E = 0;
    for (int i=0; i<R; i++){
        for (int j=0; j<C; j++){
            E += - Lattice(i, j) * Lattice(i, (j + 1) % C);
            E += - Lattice(i, j) * Lattice((i + 1) % R, j);
        };
    };
    return E;
}

void mc_cycle(arma::mat& Lattice, double& T, double& e, double& m, double& Cv, double& chi){
    /*
    Runs a single Monte-Carlo cycle of the system.

    Arguments:
        Lattice: arma::mat
            Object describing the lattice of spins.
        T: double
            The temperature of the system.
        e: double
            The energy per spin of the system.
        m: double
            The magnetisation per spin of the system.
        Cv: double
            The specific heat capacity of the system.
        chi: double
            The susceptibility of the system.
    XXX Veit ikke mer om hvordan jeg skal sette dette her opp
    */

    double E = calc_E(Lattice);
	double E_sum = 0;
	double E_2_sum = 0;
	int M = arma::accu(Lattice);
	double M_sum = 0;
	double M_2_sum = 0;
    int L = Lattice.n_rows;
    int N = L * L;

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> unifN(0, N - 1); // distribution in range [0, N - 1] for random index of attempted flip
	std::uniform_real_distribution<> flip(0.0, 1.0);						  // distribution in range [0, 1] for chance of flip

	double beta = 1 / T;
	arma::vec DEs(5, arma::fill::zeros);
	for (int i = 0; i <= 4; i += 1)
	{
		int j = -8 + i * 4;
		DEs(i) = exp(-beta * j);
	}

    for (int mc = 0; mc < N; mc++)
	{
        int idx = unifN(rng);
		int y = idx % L;
		int x = (idx - y) / L;
		int sum_neighours = Lattice((x + L - 1) % L, y) + Lattice((x + 1) % L, y) + Lattice(x, (y + L - 1) % L) + Lattice(x, (y + 1) % L);
		int DeltaE = sum_neighours * Lattice(x, y) * 2;
		int DeltaEidx = DeltaE / 4 + 2;
		double Boltzmann = DEs[DeltaEidx];

		if (Boltzmann > flip(rng))
		{
			Lattice(x, y) *= -1;
			E += DeltaE;
			M += Lattice(x, y) * 2;
		}

        E_sum += E;
		E_2_sum += E * E;
		M_sum += abs(M);
		M_2_sum += M * M;
	}

	e = E_sum / (N * N);
	E_2_sum /= N * N;
	Cv = beta / T * (E_2_sum - e * e);
	m = M_sum / (N * N);
	M_2_sum /= N * N;
	chi = beta * (M_2_sum - m * m);
}

arma::mat mc_run_culm(int L, int M, double T, std::string method="random", int burnin=0){
    /*
    Runs a number of Monte-Carlo cycles and saves the output at each cycle.

    Arguments:
        L: int
            Size of the lattice.
        M: int
            Number of Monte-Carlo cycles.
        T: double
            The temperature of the system.
        method: std::string
            Method used to initialise the lattice.
            Can be one of "random", "lowest" or "highest".
            Default is set as "random".
    Returns:
        Data: arma::mat
            Matrix containing all of the data from the Monte-Carlo cycles:
                Line 0: Values of the energy per spin given by the Monte-Carlo cycle.
                Line 1: Culmulative average of the energy per spin.
                Line 2: Values of the magnetisation per spin given by the Monte-Carlo cycle.
                Line 3: Culmulative average of the magnetisation per spin.
                Line 4: Values of the specific heat capacity given by the Monte-Carlo cycle.
                Line 5: Culmulative average of the specific heat capacity.
                Line 6: Values of the susceptibility given by the Monte-Carlo cycle.
                Line 7: Culmulative average of the susceptibility.
    */

    double e;
    double m;
    double Cv;
    double chi;
    arma::mat Data(8, M);
    arma::mat Lattice = make_sys(L, method);

    for (int i = 0; i < M; i++){
        mc_cycle(Lattice, T, e, m, Cv, chi);
        Data(0, i) = e;
        Data(2, i) = m;
        Data(4, i) = Cv;
        Data(6, i) = chi;
        if (i <= burnin){
            Data(1, i) = Data(0, i);
            Data(3, i) = Data(2, i);
            Data(5, i) = Data(4, i);
            Data(7, i) = Data(6, i);
        }
        else{
            Data(1, i) = (Data(1, i - 1) * i + Data(0, i)) / (i + 1 - burnin);
            Data(3, i) = (Data(3, i - 1) * i + Data(2, i)) / (i + 1 - burnin);
            Data(5, i) = (Data(5, i - 1) * i + Data(4, i)) / (i + 1 - burnin);
            Data(7, i) = (Data(7, i - 1) * i + Data(6, i)) / (i + 1 - burnin);
        }
    }
    return Data;
}


void mc_run(int L, int M, double T, double& e_ave, double& m_ave, double& Cv_ave, double& chi_ave, std::string method="random", int burnin=0){
    /*
    Runs a number of Monte-Carlo cycles and gives the final output.

    Arguments:
        L: int
            Size of the lattice.
        M: int
            Number of Monte-Carlo cycles.
        T: double
            The temperature of the system.
        method: std::string
            Method used to initialise the lattice.
            Can be one of "random", "lowest" or "highest".
            Default is set as "random".
        burnin: int
            Number of Monte-Carlo cycles to skip over after initialisation.

    XXX
    */

    double e;
    double m;
    double Cv;
    double chi;
    e_ave = 0;
    m_ave = 0;
    Cv_ave = 0;
    chi_ave = 0;
    arma::mat Data(8, M);
    arma::mat Lattice = make_sys(L, method);

    for (int i = 0; i < M; i++){
        mc_cycle(Lattice, T, e, m, Cv, chi);
        if (i >= burnin){
            e_ave += e;
            m_ave += m;
            Cv_ave += Cv;
            chi_ave += chi;
        }
    };
    e_ave /= (M - burnin);
    m_ave /= (M - burnin);
    Cv_ave /= (M - burnin);
    chi_ave /= (M - burnin);
}

arma::mat mc_e_prob(arma::mat& Lattice, double T, int M){
    /*
    Gives the probability density function as estimated by Monte-Carlo runs.

    Arguments:
        Lattice: arma::mat
            Object describing the lattice of spins.
        T: double
            The temperature of the system.
        M: int
            Number of Monte-Carlo cycles.
    XXX Godt spørsmål her om vi skal gjøre det på denne måten eller ha en funksjon for en syklus og legge dem alle sammen til slutt?
    */

    int L = Lattice.n_rows;
    int N = L * L;
    double E = calc_E(Lattice);
    arma::vec E_vals = arma::linspace(-2 * N, 2 * N, N + 1);
    E_vals /= N;
    arma::vec E_density(N + 1, arma::fill::zeros);

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> unifN(0, N - 1); // distribution in range [0, N - 1] for random index of attempted flip
	std::uniform_real_distribution<> flip(0.0, 1.0);						  // distribution in range [0, 1] for chance of flip

	double beta = 1 / T;
	arma::vec DEs(5, arma::fill::zeros);
	for (int i = 0; i <= 4; i += 1)
	{
		int j = -8 + i * 4;
		DEs(i) = exp(-beta * j);
	}

    for (int mc = 0; mc < N * M; mc++)
	{
        int idx = unifN(rng);
		int y = idx % L;
		int x = (idx - y) / L;
		int sum_neighours = Lattice((x + L - 1) % L, y) + Lattice((x + 1) % L, y) + Lattice(x, (y + L - 1) % L) + Lattice(x, (y + 1) % L);
		int DeltaE = sum_neighours * Lattice(x, y) * 2;
		int DeltaEidx = DeltaE / 4 + 2;
		double Boltzmann = DEs[DeltaEidx];

		if (Boltzmann > flip(rng))
		{
			Lattice(x, y) *= -1;
			E += DeltaE;
		}
        int Eidx = (E + 2 * N) / 4;
        E_density(Eidx) += 1;
	}
    E_density /= N * M;
    arma::mat E_prob;
    E_prob.insert_cols(0, E_vals);
    E_prob.insert_cols(1, E_density);
    return E_prob;
}

void multi_mc(int L, int M, int R, double T, double& e_ave, double& m_ave, double& Cv_ave, double& chi_ave, int burnin=0){
    /*
    Her skal lages en funksjon som løper over fleire initialiseringer og regner gjennomsnittet (av gjennomsnittene).
    Den skal ha en opsjon om den er parallelisert eller ikke.
    */
    double e;
    double m;
    double Cv;
    double chi;
    e_ave = 0;
    m_ave = 0;
    Cv_ave = 0;
    chi_ave = 0;
    #pragma omp parallel for
    for (int i = 0; i <= R; i++){
        mc_run(L, M, T, e, m, Cv, chi, "random", burnin);
        e_ave += e;
        m_ave += m;
        Cv_ave += Cv;
        chi_ave += chi;
    }
    e_ave /= R;
    m_ave /= R;
    Cv_ave /= R;
    chi_ave /= R;
}