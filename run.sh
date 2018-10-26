#!/bin/bash

# run.sh

rm -rf results bin
mkdir -p results bin

make clean -s
make -s

echo "generating graph"
python _data/make-graph.py

echo "running various versions"

# original c++ version
# ./bin/gtf > results/graphtv

# python version

# pure python implementation
python main.py > results/python
python reference/sfl.py \
    --validate results/python \
    --node-path _data/n  \
    --edge-path _data/e

python reference/sfl.py \
    --node-path _data/n \
    --edge-path _data/e > results/reference

# # FFA c++ version
./bin/ffa > results/ffa
python reference/sfl.py \
    --validate results/ffa \
    --node-path _data/n  \
    --edge-path _data/e


# echo "==========================="
# echo "validate gtf ffa"
# python validate.py results/graphtv results/ffa

# echo "==========================="
# echo "validate gtf reference"
# python validate.py results/graphtv results/reference

# echo "==========================="
# echo "validate ffa reference"
# python validate.py results/ffa results/reference
