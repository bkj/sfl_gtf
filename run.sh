#!/bin/bash

# run.sh

rm -rf results
mkdir -p results

make clean
make

# original c++ version
./gtf > results/gtf

# FFA c++ version
./ffa > results/ffa

python validate.py results/gtf results/ffa

# python version
python reference/sfl.py \
    --node-path _data/nodes.txt \
    --edge-path _data/edges.txt > results/reference

# check results
python validate.py results/gtf results/reference
python validate.py results/ffa results/reference
