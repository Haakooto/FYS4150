g++ main.cpp -o main.out
if [ "$3" ]; then
    ./main.out "$1" "$2" "optim"
else
    ./main.out "$1" "$2"
fi
