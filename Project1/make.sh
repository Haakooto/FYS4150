g++ main.cpp -o main.out
if [ "$2" ]; then
    ./main.out "$1" "optim"
else
    ./main.out "$1"
fi
