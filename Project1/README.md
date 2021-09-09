# Project 1 - Numerically solving 1D Poisson equation in C++

Sara Pernille Jensen, Alec Elías Sigurðarson, Håkon Olav Torvik

This repository contains our code for project 1 of Fys4150 fall 2021. This README deatils the code, including how to run it.

## Running code
### [make.sh](https://github.com/Haakooto/FYS4150/blob/main/Project1/make.sh)
Compiles and links the C++ code, as well as running it for the specifed input arguments. Subdirectory `datas` required in directory for output datafiles.
- n: int. Max log(N) for run algorithm for
- m: int. Number of repetitions for each n
- algorithm: any,  optional. If present, runs optimised Thomas algorithm. Otherwise run general Thomas algorithm.

Example: "`$ . make.sh 6 3 o`" runs optimised algorithm 3 times for each $n \in [1, 6]$ .

### [plot.py](https://github.com/Haakooto/FYS4150/blob/main/Project1/plot.py)
Makes various plots of the data, dependant on input arguments. Full description of each plotting function is given below.
- n: int. Same as for make.sh
- function: string, optional. Name of plotting function.  Defaults to 'num_sol'.
- algorithm: any, optional. Same as for make.sh. Not used in all plotting functions.

Example: "`$ python plot.py 6 timing`"

## Code Overview
File-by-file rundown and explanation of each function present in code.
### [main.cpp](https://github.com/Haakooto/FYS4150/blob/main/Project1/main.cpp)
Starting point for every run. Compiled and run from make.sh. Loops over n and m to run specified algorithm.

### [thomas.hpp](https://github.com/Haakooto/FYS4150/blob/main/Project1/thomas.hpp)
~~Technically wrong use of header files, but it works~~.

General Thomas algorithm. Called by main.cpp.

#### Functions:
- Thomas(): Allocates required memory to arrays, has them initialized, runs algorithm, writes output to file, frees memory.
- Fwd_Bkwd_sub(): Performs algorithm.
- init_arrays(): Initialises arrays

### [optimized_thomas.hpp](https://github.com/Haakooto/FYS4150/blob/main/Project1/optimized_thomas.hpp)
Exactly same as thomas.hpp, but for the optimized algorithm. Is faster and requires less memory.

### [utils.hpp](https://github.com/Haakooto/FYS4150/blob/main/Project1/utils.hpp)
Utility functions.
#### Functions:
- f(): given source function.
- analytic_sol(): analytical solution
- write_full(): writes entire x-domain, analytic solution, numeric solution, absolute error and relative error at every point $x_i$ to file.
- write_limited(): writes n (log10(N)), time spent on algorithm, and maximum relative error to file for all n and m.
- rel_err(): calculates relative error at every point, and return maximum relative error
- abs_err(): calculates absolute error at every point.

### [plot.py](https://github.com/Haakooto/FYS4150/blob/main/Project1/plot.py)
Using library Plotly to make very nice graphs.
#### Functions():
- f() and data(): Loads dataframes. Do not attempt to understand how, what, why nor when.
- num_sol(): Plots numerical and analytical solution for a single n
- first_few_num_sol(): Same as num_sol() for for multiple n in same plot.
- absolute_error(): Plots absolute error at every point for single n.
- relative_error(): Plots relative error at every point for single n.
- max_error(): Plots log10(MRE) for multiple n together with a fitted line.
- both_max_error(): Same as max_error(), but for both general and optimized algorithm in same plot.
- timing(): Plots log10(time spent on algorithm) for both general and optimized algorithm. If multiple runs for each n, calculate mean and variance, plot mean with error bars. Also a line fitted to the data.