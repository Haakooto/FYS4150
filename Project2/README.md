# Project 2 - Implementing Jacobi's Rotation Method for solving an eigenvalue problem in C++
Sara Pernille Jensen, Alec Elías Sigurðarson, Håkon Olav Torvik

This repository contains our code for Project 2 of FYS4150 fall 2021.

## Running code
### make.sh
Write ". make.sh" in terminal to run. 
Compiles and links the C++ code, and runs the main function from main.cpp. Runs all unit tests, as well as writing all the required data to plot the graphs. 
The python file must be run separately to generate the graphs. 
Subdirectory "results" required in directory to save datafiles in. 

### plot.py
Write "plot.py" in terminal to generate graphs. 
load: takes the correct data-file saved in "results".
Change the value of N manually in the script to look at results for different values. 

#### Functions:
iterations(): Plots the number of transformations as a function of N. 
solutions(): Plots the three first solutions obtained from the Jacobi method and the analytic solutions. 
main(): XXX

## Code Overview
Explanations of each file, including its funstions.

### main.cpp
Needed to solve the equations and write the files. Compiled and run from make.sh. 

#### Functions:
max_offdiag_symmetric(): Finds the value of position of the largest off-diagonal element of the matrix.
Rotation(): Performs the rotation and updates the values of the elements of 
Jacobi(): Calls on the two functions above until the largest off-diagonal element has a value below the defined tolerance. 
test_max_offdiag(): Unit test to check that max_offdiag_symmetric() returns the correct value and position. 
test_analyticity(): Unit test to check that the analytic solutions correspond to armadillo's solutions to the same problem. 
test_algorithm(): Unit test to check that the solution obtained from Jacobi's rotation method corresponds with the analytic solution. 
tests(): Runs the three abovementioned unit tests. 
run_Jacobi(): Runs the Jacobi rotation method for a given matrix (using Jacobi()), sorts the solution and writes this to a file. 
main(): Runs the Jacobi algorithm for specified values of N and maximum number of iterations. 

### utils.cpp
Utility funcitons.

#### Functions:
make_A(): Constructs the matrix A for a given dimensionality and h. 
analytic_eigval(): Calculates the eigenvalues analytically.  
analytic_eigvec(): Calculates the eigenvectors analytically and returns the normalised vectors. 
analytic_solutions(): Generates the complete solution by calling on both analytic_eigval() and analytic_eigvec(). 
max_diff(mat): Calculates the maximum difference between the elements in two matrices. Used for unit testing. 
max_diff(vec): Calculates the maximum difference between the elements in two vectors. Used for unit testing.
sort_mat_by_vec(): Sorts the columns of the matrix in ascending order of their corresponding eigenvalues. 
write_to_file(): Creates a new file and writes out the eigenvalues and eigenvectors provided by the solution. The boundary conditions are added manually. 


### utils.hpp
Includes the functions needed. 
XXX


