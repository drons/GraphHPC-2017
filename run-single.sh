#!/bin/sh
GSIZE=10
#export OMP_NUM_THREADS=56
#export KMP_AFFINITY=scatter,granularity=thread

./gen_RMAT -s $GSIZE &&
./solution -nIters 2 -in rmat-$GSIZE &&
./gen_valid_info -in rmat-$GSIZE &&
./validation -ans rmat-$GSIZE.ans -res rmat-$GSIZE.res &&

./gen_random -s $GSIZE &&
./solution -nIters 2 -in random-$GSIZE &&
./gen_valid_info -in random-$GSIZE &&
./validation -ans random-$GSIZE.ans -res random-$GSIZE.res &&
exit 0
