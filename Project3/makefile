compiletest:
	g++ --std=c++20 test.cpp -larmadillo src/* -I include/ -o test.exe

all:
	compiletest

wafti:
	g++ --std=c++20 wafti.cpp -larmadillo -o test.exe
	./test.exe