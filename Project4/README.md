# Project4 - The Ising model
Sara Pernille Jensen, Alec Elías Sigurðarson, Håkon Olav Torvik

This repository contains our code for Project 4 of FYS4150 fall 2021.

# Ising Model
We have studied the ising model using Monte-Carlo Markov chains to sample spin configurations through the Metropolis algorithm.

# Code 
The central algorithm and support functions for this is found in the file `ìsing_model.cpp`. Here are also several wrapper functions to run with different configurations. The various other experiment `.cpp`-files call upon one of these wrapper functions. Each file either saves a file in the `data`-directory, or prints output to terminal. The python-file `plot.py` is used to analyse and/or visualize the data produced by these. It contains a function for each experiemnt, and makes its own data where nessisary by calling on the needed experment. Some of the functions take certain commandline arguments. 

## Compiling
All `.cpp`-files must be compiled with the linear algebra library Armadillo. Several functions are parallelized using OpenMP, so it is also benefitial to add the flag `-fopenmp` when compiling. The makefile includes commands for the compiling we have used. These can all be run by `$ make compile_all`. 

## Running code
The makefile also includes the commands to run the `plot.py` to generate the results for each problem. For instance, `$ make ex_4` will print the tables of analytical values to the terminal, latex-ready. `$ make all` will generate all the plots used in the rapport, and print the numerical results where appropriate. 

Alternatively, `plot.py` can be run directly from command line. Then the name of the desired function must be passed as a commandline argument, as well as the possible other arguments the function takes. These can be found the the respective docstrings.  

Listed below are all the C++ experiment files, the corrosponding python function to analyse/visualize the data, and a short description of what it does.
- `analytic.cpp` with `analytic()` Makes table of numerical vs analytical values and prints to terminal
- `burning.cpp` with `burntime()`  Makes plot of energy and magnetization as function of MC-cycles
- `pdf.cpp` with `pdf()` Makes histogram of energy distribution
- `MC_Hammer.cpp` with `paralympics()`  Prints time taken by algorithm with and without optimizationflag and parallelization
- `tempting.cpp` with `run_temps()`  Make plot of energy, magnetization, heat capacity and susceptibility as function of temperature for different lattice sizes. `critical_temp()` can further be used to fit the heat capacity and susceptibility and from this determine the critical temperature for an infinite lattice size. 

## Dependencies
As noted already, Armadillo and OpenMP is used for the C++ code. The following modules not always in the standard python installation are used by `plot.py`:
- numpy
- scipy
- pandas
- plotly
- subprocess
- uncertainties
