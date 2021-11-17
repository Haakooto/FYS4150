compile:
	g++ -std=c++11 test.cpp -larmadillo -o test.out

run:
	./test.out

all: compile run

fire:
	g++ burnin.cpp -o burn -larmadillo

ex_5:
	python3 plot.py burntime burn_1 1 100 True
	python3 plot.py burntime burn_24 2.4 300 True

compile_analytic:
	g++ -std=c++11 analytic.cpp -larmadillo -o analytic.out

run_analytic:
	./analytic.out

analytic_sol: compile_analytic run_analytic
