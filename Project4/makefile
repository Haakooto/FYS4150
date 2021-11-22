analytic:
	g++ -std=c++11 analytic.cpp -larmadillo -o analytic.out

paralytic_mac:
	g++ -std=c++11 analytic.cpp -Xpreprocessor -fopenmp -larmadillo -o paralytic.out -lomp

paralytic_ubuntu:
	g++ -std=c++11 analytic.cpp -fopenmp -larmadillo -o paralytic.out

fire:
	g++ -std=c++11 burning.cpp -o burn.out -larmadillo

pdf:
	g++ -std=c++11 pdf.cpp -larmadillo -o pdf.out

paradef:
	g++ pdf.cpp -o pdf.out -fopenmp -larmadillo

cant_touch_this:
	g++ MC_Hammer.cpp -o time.out -larmadillo -fopenmp
	g++ MC_Hammer.cpp -o optime.out -larmadillo -fopenmp -O3

tempting:
	g++ tempting.cpp -o tempting.out -larmadillo -fopenmp -O3

compile_all: analytic paralytic_ubuntu fire pdf paradef cant_touch_this tempting

ex_4:
	python3 plot.py analytic 1 23456 100 True
	python3 plot.py analytic 2.4 23456 100 True

ex_5:
	python3 plot.py burntime burn 1 1000
	python3 plot.py burntime burn 2.4 3000

ex_6:
	./pdf.out 10000 1
	python3 plot.py pdf

ex_7:
	python3 plot.py paralympics

ex_8:
	python3 plot.py run_temps temps 2.1 2.5 30 1000000 True '40,60,80,100'
	python3 plot.py run_temps temps_narrow 2.25 2.325 15 5000000 True '40,60,80,100'

ex_9:
	python3 plot.py critical_temp temps_narrow 2.25 2.325 15 '40,60,80,100'

all:
	ex_4 ex_5 ex_6 ex_7 ex_8 ex_9 