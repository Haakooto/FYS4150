#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "utils.hpp"
#include "thomas.hpp"
#include "optimized_thomas.hpp"

using namespace std;

int main(int argc, char *argv[]){
    int n;
    n = atoi(argv[1]);

    // To not wait forever for writing full output, cap such writing for low n
    int cap = 4;

    void (*Method)(int, int, double*, double*) = NULL;
    string name;
    if (argc == 3){
        name = "_optim_";
        Method = Optimized_Thomas;
    } else {
        name = "";
        Method = Thomas;
    }
    int ns[n];
    double mres[n];
    double times[n];

    for (int i = 1; i < n + 1; i++){
        cout << "starting " + name + "Thomas algorithm for n = " << i << endl;
        Method(i, cap, times, mres);
        ns[i - 1] = i;
    }
    write_limited(n, ns, times, mres, name);
    return 0;
}
