#!/bin/bash

rm *.o
rm run_capelec

make
# export OMP_NUM_THREADS=4

./run_capelec $1 0
