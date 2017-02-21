#!/bin/sh
GSIZE=17
export OMP_NUM_THREADS=56
export KMP_AFFINITY=scatter,granularity=thread
./solution -nIters 1 -in random-$GSIZE
