#include <iostream>
#include <string>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <armadillo>
#include <random>
#include <omp.h>


arma::vec make_de(double beta){
    /*
    Arguments:
        beta: double
            1 over temperature of system

    Returns
        DEs: arma::vec
            Vector with Boltzmann probability factors for allowed energy changes.
    */
    arma::vec DeltaE(5, arma::fill::zeros);
    for (int i = 0; i <= 4; i += 1)
    {
        int j = -8 + i * 4;
        DeltaE(i) = exp(-beta * j);
    }
    return DeltaE;
}

arma::mat make_sys(int L, std::string method="random"){
    /*
    Initialises an Ising model lattice.

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
        arma::arma_rng::set_seed_random();
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
    Calculates the energy of a system.

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

void mc_cycle(arma::mat& Lattice, arma::vec& DEs, double& E_sum, double& M_sum, double& E_sq, double& M_sq, int seed){
    /*
    Runs a single Monte-Carlo cycle of the system.
    Core of the matter

    Arguments:
        Lattice: arma::mat
            Object describing the lattice of spins.
        DEs: arma::vec
            Vector with Boltzmann probability factors for allowed energy changes.
        E_sum: double
            The energy per spin of the system.
        M_sum: double
            The magnetisation per spin of the system.
        E_sq: double
            The specific heat capacity of the system.
        M_sq: double
            The susceptibility of the system.
        seed: int
            seed for random number generator.
            MUST NEVER BE THE SAME SEED FOR A RUN.
            For example, use loop indexer when calling mc_cycle in a loop
    */

    double E = calc_E(Lattice);
	E_sum = 0;
	E_sq = 0;
	int M = arma::accu(Lattice);
	M_sum = 0;
	M_sq = 0;
    int L = Lattice.n_rows;
    int N = L * L;

	std::mt19937 rng(seed);
	std::uniform_int_distribution<std::mt19937::result_type> unifN(0, N - 1); // distribution in range [0, N - 1] for random index of attempted flip
	std::uniform_real_distribution<> flip(0.0, 1.0);						  // distribution in range [0, 1] for chance of flip

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
		E_sq += E * E;
		M_sum += abs(M);
		M_sq += M * M;
	}


	E_sum /= N;
	E_sq /= N;

	M_sum /= N;
	M_sq /= N ;
}

arma::mat mc_run_cuml(int L, int M, double T, std::string method="random", int burnin=0){
    /*
    Runs a number of Monte-Carlo cycles and saves the output at each cycle.
    Used for calculating burntime

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
                Line 1: Cumulative average of the energy per spin.
                Line 2: Values of the magnetisation per spin given by the Monte-Carlo cycle.
                Line 3: Cumulative average of the magnetisation per spin.
    */

    double e;
    double m;
    double EE;
    double MM;
    double N = L * L;
    arma::mat Data(4, M);
    arma::mat Lattice = make_sys(L, method);
    arma::vec DEs = make_de(1 / T);


    for (int i = 0; i < M; i++){
        mc_cycle(Lattice, DEs, e, m, EE, MM, i);
        Data(0, i) = e / N;
        Data(2, i) = m / N;
        if (i <= burnin){
            Data(1, i) = Data(0, i);
            Data(3, i) = Data(2, i);
        }
        else{
            Data(1, i) = (Data(1, i - 1) * i + Data(0, i)) / (i + 1 - burnin);
            Data(3, i) = (Data(3, i - 1) * i + Data(2, i)) / (i + 1 - burnin);
        }
    }
    return Data;
}

void mc_run(int L, int M, arma::vec DEs, double& e_ave, double& m_ave, double& Cv_ave, double& chi_ave, std::string method="random", int burnin=0){
    /*
    Runs a number of Monte-Carlo cycles and gives the final output.

    Arguments:
        L: int
            Size of the lattice.
        M: int
            Number of Monte-Carlo cycles.
        DEs: arma::vec
            Vector with Boltzmann probability factors for allowed energy changes.
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
    double N = L * L;
    e_ave = 0;
    m_ave = 0;
    Cv_ave = 0;
    chi_ave = 0;
    arma::mat Lattice = make_sys(L, method);
    for (int mc = 0; mc < M; mc++){
        mc_cycle(Lattice, DEs, e, m, Cv, chi, mc);
        if (mc >= burnin){
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
    Cv_ave = 1 / N * (Cv_ave - pow(e_ave, 2));
    chi_ave = 1 / N * (chi_ave - pow(m_ave, 2));
    e_ave /= N;
    m_ave /= N;
}

arma::mat mc_e_prob(arma::mat& Lattice, double T, int M, int burnin=0){
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

	arma::vec DEs = make_de(1 / T);

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
        if (mc >= burnin){
            int Eidx = (E + 2 * N) / 4;
            E_density(Eidx) += 1;
        }
	}
    E_density /= N * M;
    arma::mat E_prob;
    E_prob.insert_cols(0, E_vals);     //the possible e-values
    E_prob.insert_cols(1, E_density);   // the likelyhood of each e-value
    return E_prob;
}

void multi_mc(int L, int M, int R, double T, arma::vec& data, std::string method="random", int burnin=0, bool para=false)
{
     /*
    Run R cocurrent mcmc with the same temperature, return mean of calculated quantities

    Arguments:
        L: int
            Size if lattice
        M: int
            Number of MCMC cycles
        R: int
            Number of cocurrent runs for same temperature
        T: double
            Temperature of system
        data: arma::vec
            data-containing vector
                index 0: energy / spin
                index 1: magnetizaton / spin
                index 2: heatcapacity
                index 3: suceptibility
                index 4: std of index 0 from R runs
                index 5: std of index 1 from R runs
                index 6: std of index 2 from R runs
                index 7: std of index 3 from R runs
        method: string
            Initialization method
        burnin: int
            Number of burn-in cycles
    XXX
    */

    arma::vec e_vec(R, arma::fill::zeros);
    arma::vec m_vec(R, arma::fill::zeros);
    arma::vec Cv_vec(R, arma::fill::zeros);
    arma::vec chi_vec(R, arma::fill::zeros);

    arma::vec DEs = make_de(1 / T);
    if (para){
        #pragma omp parallel for
            for (int i = 0; i < R; i++)
            {
                double e, m, Cv, chi;
                mc_run(L, M, DEs, e, m, Cv, chi, method, burnin);
                e_vec(i) = e;
                m_vec(i) = m;
                Cv_vec(i) = Cv;
                chi_vec(i) = chi;
            }
    } else {
            for (int i = 0; i < R; i++)
            {
                double e, m, Cv, chi;
                mc_run(L, M, DEs, e, m, Cv, chi, method, burnin);
                e_vec(i) = e;
                m_vec(i) = m;
                Cv_vec(i) = Cv;
                chi_vec(i) = chi;
            }
    }

    Cv_vec /= (T * T);
    chi_vec /= T;

    data(0) = arma::mean(e_vec);
    data(1) = arma::mean(m_vec);
    data(2) = arma::mean(Cv_vec);
    data(3) = arma::mean(chi_vec);

    data(4) = arma::stddev(e_vec)/sqrt(R);
    data(5) = arma::stddev(m_vec)/sqrt(R);
    data(6) = arma::stddev(Cv_vec)/sqrt(R);
    data(7) = arma::stddev(chi_vec)/sqrt(R);
}

arma::mat multi_prob(int L, int M, int R, double T, int burnin=0){
    /*
    Runs mc_e_prob R times and averages result. Is parallellized over R
    Arguments:
        L: int
            Think you know this shit by now
        M: int
            Like, seriously, there is nothing new under the LED-bubls
        R: int
            Maybe I should take up poetry?
        T: double
            Would keep me busy for a while
        burnin: int
            Actually, that is a great idea
    */
    arma::mat Total(L * L + 1, 2, arma::fill::zeros);
    arma::mat sum;
    #pragma omp parallel for

        for (int i = 0; i < R; i++){
            arma::mat Lattice = make_sys(L, "random");
            sum = mc_e_prob(Lattice, T, M, burnin);
            Total += sum;
        }

    Total /= R;
    return Total;
}

arma::mat thepoem(int L, int M, int R, double T, int burnin=0){
    int N = L * L;
    // What do I do? Not even Enya can answer that one...
    arma::mat Total(N + 1, 2, arma::fill::zeros);
    // Stay beautiful, keep it ugly
    arma::mat sum;
    // As with anything creative, change is inevitable
    for (int i = 0; i < R; i++){
        // Think I should exist in another file...
        arma::mat Lattice = make_sys(L, "random");
        // Would you destroy something perfect in order to make it beautiful
        sum = mc_e_prob(Lattice, T, M, burnin);
        // Perhaps I am already dead
        Total += sum;
        // At the very least living in another scope
    }
    // I wanted to become God, destroyer of death
    Total /= R;
    // I at least attempted to give back the sum of it all
    return Total;
    // Now, I walk gentle into that good night
}
