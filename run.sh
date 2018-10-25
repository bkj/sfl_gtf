#!/bin/bash

# run.sh

rm -rf results
mkdir -p results

make clean -s
make -s

# echo "generating graph"
# python _data/make-graph.py
python num_cc.py

# --
# BKJ c++ version

echo "============================================="
echo "c++ (bkj)"

./bkj _data/nodes.txt _data/edges.txt > results/bkj

python reference/sfl.py \
    --validate results/bkj \
    --node-path _data/nodes.txt  \
    --edge-path _data/edges.txt

# --
# FFA C++ version

echo "============================================="
echo "c++ (ffa)"

./ffa _data/nodes.txt _data/edges.txt > results/ffa

python reference/sfl.py \
    --validate results/ffa \
    --node-path _data/nodes.txt  \
    --edge-path _data/edges.txt

# --
# python version

echo "============================================="
echo "python"

python reference/sfl.py \
    --node-path _data/nodes.txt \
    --edge-path _data/edges.txt > results/reference