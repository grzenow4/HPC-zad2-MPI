#!/bin/bash

if [ "$1" = "example" ]; then
    mpiexec -n 16 ./matmul -a ../tests/sparse_matrix_a -b ../tests/sparse_matrix_b -t 2D -v -g 20
elif [ "$1" = "paper" ]; then
    mpiexec -n 12 ./matmul -a ../tests/paper_a_b -b ../tests/paper_a_b -t 3D -l 3 -v -g 2
else
    echo "Usage: $0 [example|paper]"
    exit 1
fi
