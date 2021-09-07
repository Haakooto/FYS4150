n="$1"
if [ "$2" ]; then
    method="$2"
else
    method="Thomas"
fi

g++ test.cpp -o main.out
./main.out "$n" "$method"
python3 plot.py "$n"
