compile:
	g++ -std=c++11 test.cpp -larmadillo -o test.out

run:
	./test.out

all: compile run

fire:
	g++ burnin.cpp -o burn -larmadillo

compile_analytic:
	g++ -std=c++11 analytic.cpp -larmadillo -o analytic.out

run_analytic:
	./analytic.out

analytic_sol: compile_analytic run_analytic
