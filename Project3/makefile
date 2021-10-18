compiletest:
	g++ --std=c++20 test.cpp -larmadillo src/* -I include/ -o test.exe

all:
	compiletest

wafti:
	g++ --std=c++20 wafti_old.cpp -larmadillo -o test.exe
	./test.exe

main:
	g++ -std=c++20 main.cpp src/* -I include/ -o main.out -larmadillo
	# ./main.out

run:
	./main.out
