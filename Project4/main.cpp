#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <random>

using namespace std;

int main()
{

	int L = 2;
	int N = pow(L, 2);
	double T = 1.;

	std::random_device dev;
	std::mt19937 rng(dev());
	std::uniform_int_distribution<std::mt19937::result_type> unifN(0, N - 1); // distribution in range [0, 1]
	std::uniform_real_distribution<> flip(0.0, 1.0);						  // distribution in range [0, 1]

	double beta = 1 / T;
	arma::vec DEs(5, arma::fill::zeros);
	for (int i = 0; i <= 4; i += 1)
	{
		int j = -8 + i * 4;
		DEs(i) = exp(-beta * j);
	}

	int M = 500;

	double E = -8;
	double E_sum = E;
	double E_2_sum = E * E;
	int m = 4;
	double m_sum = m;
	double m_2_sum = m * m;

	arma::mat Lattice(L, L, arma::fill::ones);
	// Lattice[1, 0] = -1;

	for (int mc = 0; mc < M * N; mc++)
	{
		int idx = unifN(rng);
		int y = idx % L;
		int x = (idx - y) / L;
		// cout << (x + 1) % L << " " << (x - 1) % L << " " <<(y - 1) % L << " " <<(y + 1) % L << endl;
		int sum_neighours = Lattice((x + L - 1) % L, y) + Lattice((x + 1) % L, y) + Lattice(x, (y + L - 1) % L) + Lattice(x, (y + 1) % L);
		int DeltaE = sum_neighours * Lattice(x, y) * 2;
		int DeltaEidx = DeltaE / 4 + 2;
		double Boltzmann = DEs[DeltaEidx];
		// DeltaE /= 2 + 2;
		// DeltaE += 2;
		// cout << sum_neighours << " " << DeltaEidx << " " << Boltzmann << " " << DeltaE << endl;
		// cout << idx << " " << x << " " << y << endl;

		if (Boltzmann > flip(rng))
		{
			Lattice(x, y) *= -1;
			E += DeltaE;
			m += Lattice(x, y) * 2;
		}
		E_sum += E;
		E_2_sum += E * E;
		m_sum += abs(m);
		m_2_sum += m * m;
	}

	double E_avg = E_sum / (M * N * N);
	E_2_sum /= M * N * N;
	double Cv = beta / T * (E_2_sum - E_avg * E_avg);
	double m_avg = m_sum / (M * N * N);
	m_2_sum /= M * N * N;
	double chi = beta * (m_2_sum - m_avg * m_avg);

	cout << E_avg << " " << Cv << " " << m_avg << " " << chi << endl;

	return 0;
}