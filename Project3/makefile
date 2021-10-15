compiletest:
	g++ --std=c++20 test.cpp -larmadillo src/* -I include/ -o test.exe

all:
	compiletest

wafti:
	g++ --std=c++20 wafti.cpp -larmadillo -o test.exe
	./test.exe

main:
	g++ -std=c++20 main.cpp -larmadillo src/* -I include/ -o main.out
	./main.out

run:
	./main.out
