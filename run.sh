#!/bin/bash

# run.sh

rm -rf results
mkdir -p results

make clean -s
make -s

echo "generating graph"
python _data/make-graph.py

# original c++ version
./gtf > results/gtf

# FFA c++ version
./ffa > results/ffa

# python version
python reference/sfl.py \
    --node-path _data/nodes.txt \
    --edge-path _data/edges.txt > results/reference

echo "==========================="
echo "validate gtf ffa"
python validate.py results/gtf results/ffa

echo "validate gtf reference"
python validate.py results/gtf results/reference

echo "validate ffa reference"
python validate.py results/ffa results/reference
