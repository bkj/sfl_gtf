#!/bin/bash

# run.sh

rm -rf results
mkdir -p results

make clean
make

# c++ version
./main > results/cpp

# python version
python reference/sfl.py \
    --node-path _data/nodes.txt \
    --edge-path _data/edges.txt > results/reference

# check results
python validate.py results/cpp results/reference
