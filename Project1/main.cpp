#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "utils.hpp"
#include "thomas.hpp"

using namespace std;

int main(int argc, char *argv[]){
    int n;

    n = atoi(argv[1]);
    for (int i = 1; i < n + 1; i++){
        cout << "starting Thomas algorithm for n = " << i << endl;
        Thomas(i);
    }

    return 0;
}
