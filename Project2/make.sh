g++ main.cpp utils.cpp -o main.out -O2 -std=c++20 -larmadillo
# [ "$?" -eq 0 ] && [ -z "$3" ] && ./main.out
[ "$?" -eq 0 ] && ./main.out
