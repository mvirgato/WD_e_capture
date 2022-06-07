#!/bin/bash

make clean
make
# export OMP_NUM_THREADS=4

./bin/run_capelec $1 0
