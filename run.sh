#!/bin/bash

# run.sh

rm -rf results
mkdir -p results

make clean -s
make -s

echo "generating graph"
python _data/make-graph.py

# original c++ version
# ./gtf > results/gtf

# FFA c++ version
# ./ffa > results/ffa

# BKJ c++ versoin
./bkj > results/bkj

# python version
python reference/sfl.py \
    --node-path _data/nodes.txt \
    --edge-path _data/edges.txt > results/reference

echo "==========================="
# echo "validate gtf ffa"
# python validate.py results/gtf results/ffa

# echo "validate gtf reference"
# python validate.py results/gtf results/reference

echo "validate bkj reference"
python validate.py results/bkj results/reference

# echo "validate bkj reference"
# python validate.py results/ffa results/bkj


# # >>

# python reference/sfl.py \
#     --node-path _data/taxi_n \
#     --edge-path _data/taxi_e > results/taxi_reference

# cp _data/taxi_n _data/nodes.txt
# cp _data/taxi_e _data/edges.txt

# ./ffa > results/ffa
# ./bkj > results/bkj

# # <<