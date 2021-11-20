
analytic:
	g++ -std=c++11 analytic.cpp -larmadillo -o analytic.out

paralytic_mac:
	g++ -std=c++11 analytic.cpp -Xpreprocessor -fopenmp -larmadillo -o paralytic.out -lomp

paralytic_ubuntu:
	g++ -std=c++11 analytic.cpp -fopenmp -larmadillo -o paralytic.out

fire:	
	g++ burning.cpp -o burn.out -larmadillo

pdf:
	g++ -std=c++11 pdf.cpp -larmadillo -o pdf.out

paradef:
	g++ pdf.cpp -o pdf.out -fopenmp -larmadillo

compile_all: analytic paralytic_ubuntu fire pdf paradef

ex_5:
	python3 plot.py burntime burn 1 100 True
	python3 plot.py burntime burn 2.4 300 True

ex_4:
	python3 plot.py analytic 1 23456 100 True
	python3 plot.py analytic 2.4 23456 100 True

my_day:
	g++ tempting.cpp -fopenmp -larmadillo -o tempting.out
	python3 plot.py run_temps temp_broad_search_test1 2.1 2.4 20 200 100 True '40,60,80,100'
	echo -ne '\007'
	sleep 0.2
	echo -ne '\007'
	sleep 0.2
	echo -ne '\007'
	sleep 0.2
	echo -ne '\007'

it_happen:
	g++ tempting.cpp -fopenmp -larmadillo -o tempting.out
	python3 plot.py run_temps temp_broad_search 2.1 2.4 20 10000 100 True '40,60,80,100'
	echo -ne '\007'
	sleep 0.2
	echo -ne '\007'
	sleep 0.2
	echo -ne '\007'
	sleep 0.2
	echo -ne '\007'

plots:
	python3 plot.py run_temps temp_broad_search 2.1 2.4 20 10000 100