#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "utils.hpp"
#include "thomas.hpp"
#include "optimized_thomas.hpp"

using namespace std;

int main(int argc, char *argv[]){
    int n;  // largest log(N) to run algorithm for
    int m;  // number of repetitions for each n
    n = atoi(argv[1]);
    m = atoi(argv[2]);

    // To not wait forever for writing full output, cap such writing for low n
    int cap = 4;

    // Choosing algorithm from command line
    void (*Method)(int, int, int, int, double*, double*) = NULL;
    string name;
    if (argc == 4){
        name = "_optim_";
        Method = Optimized_Thomas;
    } else {
        name = "";
        Method = Thomas;
    }

    int ns[m * n];
    double mres[m * n];  // mre = max relative error
    double times[m * n];

    for (int i = 1; i < n + 1; i++){
        for (int j = 0; j < m; j++){
            cout << "starting " + name + "Thomas algorithm for n = " << i << endl;
            Method(i, j, m, cap, times, mres);
            ns[(i - 1) * m + j] = i;
        }
    }
    write_limited(n, m, ns, times, mres, name);
    return 0;
}
