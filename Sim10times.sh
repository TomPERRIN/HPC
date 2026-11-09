
#!/bin/bash

./simPI 19937 > PI1 &
./simPI 52149 > PI2 &
./simPI 63113 > PI3 &
./simPI 86226 > PI4 &
./simPI 12345 > PI5 &
./simPI 45678 > PI6 &
./simPI 98765 > PI7 &
./simPI 24680 > PI8 &
./simPI 13579 > PI9 &
./simPI 54321 > PI10 &

wait
