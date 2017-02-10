#!/bin/sh
GSIZE=17
export OMP_NUM_THREADS=56
export KMP_AFFINITY=scatter,granularity=thread
numactl -H
./gen_RMAT -s $GSIZE && ./solution -nIters 2 -in rmat-$GSIZE
#&& ./validation -in rmat-$GSIZE -res rmat-$GSIZE.res && echo "DONE!"
