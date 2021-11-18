
analytic:
	g++ -std=c++11 analytic.cpp -larmadillo -o analytic.out

paralytic_mac:
	g++ -std=c++11 analytic.cpp -Xpreprocessor -fopenmp -larmadillo -o paralytic.out -lomp

paralytic_ubuntu:
	g++ -std=c++11 analytic.cpp -fopenmp -larmadillo -o paralytic.out

fire:
	g++ burning.cpp -o burn.out -larmadillo

pdf:
	g++ pdf.cpp -o pdf.out -larmadillo

paradef:
	g++ pdf.cpp -o pdf.out -fopenmp -larmadillo

compile_all: analytic paralytic_ubuntu fire pdf paradef

ex_5:
	python3 plot.py burntime burn 1 100 True
	python3 plot.py burntime burn 2.4 300 True

ex_4:
	python3 plot.py analytic 1 2345 100 True
	python3 plot.py analytic 2.4 2345 100 True
# kjør manuelt, eller gjennom python
# run_analytic:
# 	./analytic.out

# analytic_sol: compile_analytic run_analytic
